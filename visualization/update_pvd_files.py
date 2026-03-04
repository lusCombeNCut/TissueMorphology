#!/usr/bin/env python3
"""
Update PVD files to reference all generated VTU/VTI files.
"""
import os
import re
from pathlib import Path

def extract_timestep(filename):
    """Extract timestep number from filename (e.g., 'results_100.vtu' -> 100)"""
    match = re.search(r'_(\d+)\.vt[ui]$', filename)
    if match:
        return int(match.group(1))
    return None

def create_pvd_content(file_list, file_type='vtu'):
    """Generate PVD XML content for a list of files."""
    # Sort files by timestep
    files_with_timesteps = []
    for f in file_list:
        timestep = extract_timestep(f)
        if timestep is not None:
            files_with_timesteps.append((timestep, f))
    
    files_with_timesteps.sort(key=lambda x: x[0])
    
    # Generate XML
    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">',
        '    <Collection>'
    ]
    
    for timestep, filename in files_with_timesteps:
        lines.append(f'        <DataSet timestep="{timestep}" group="" part="0" file="{filename}"/>')
    
    lines.append('    </Collection>')
    lines.append('</VTKFile>')
    
    return '\n'.join(lines) + '\n'

def update_pvd_files(base_dir):
    """Update all PVD files in the simulation directory."""
    base_path = Path(base_dir)
    
    # Process results_from_time_5
    time5_dir = base_path / 'results_from_time_5'
    
    if time5_dir.exists():
        print(f"\nProcessing {time5_dir}...")
        
        # Get all VTU files (excluding _tri_ files)
        vtu_files = [f.name for f in time5_dir.glob('results_*.vtu') 
                     if '_tri_' not in f.name]
        print(f"  Found {len(vtu_files)} VTU files")
        
        # Get all ECM VTI files
        vti_files = [f.name for f in time5_dir.glob('ecm3d_*.vti')]
        print(f"  Found {len(vti_files)} VTI files")
        
        # Create results.pvd
        if vtu_files:
            results_pvd = time5_dir / 'results.pvd'
            content = create_pvd_content(vtu_files, 'vtu')
            with open(results_pvd, 'w') as f:
                f.write(content)
            print(f"  ✓ Created {results_pvd} with {len(vtu_files)} entries")
        
        # Create ecm3d_results.pvd
        if vti_files:
            ecm_pvd = time5_dir / 'ecm3d_results.pvd'
            content = create_pvd_content(vti_files, 'vti')
            with open(ecm_pvd, 'w') as f:
                f.write(content)
            print(f"  ✓ Created {ecm_pvd} with {len(vti_files)} entries")
    
    # Process results_from_time_0 (check if needs update)
    time0_dir = base_path / 'results_from_time_0'
    
    if time0_dir.exists():
        print(f"\nChecking {time0_dir}...")
        
        vtu_files = [f.name for f in time0_dir.glob('results_*.vtu') 
                     if '_tri_' not in f.name]
        print(f"  Found {len(vtu_files)} VTU files")
        
        results_pvd = time0_dir / 'results.pvd'
        
        # Check if PVD exists and has correct number of entries
        if results_pvd.exists():
            with open(results_pvd, 'r') as f:
                content = f.read()
                entry_count = content.count('<DataSet')
            
            if entry_count == len(vtu_files):
                print(f"  ✓ {results_pvd} is up to date ({entry_count} entries)")
            else:
                print(f"  ⚠ {results_pvd} has {entry_count} entries but found {len(vtu_files)} VTU files")
                # Uncomment to update:
                # content = create_pvd_content(vtu_files, 'vtu')
                # with open(results_pvd, 'w') as f:
                #     f.write(content)
                # print(f"  ✓ Updated {results_pvd}")
        else:
            content = create_pvd_content(vtu_files, 'vtu')
            with open(results_pvd, 'w') as f:
                f.write(content)
            print(f"  ✓ Created {results_pvd} with {len(vtu_files)} entries")

if __name__ == '__main__':
    simulation_dir = r'\\wsl.localhost\Ubuntu-20.04\home\orlando\Thesis\testoutput\VertexCryptOrganoid\growth_20260210_102802'
    
    print("=" * 60)
    print("Updating PVD files for vertex organoid simulation")
    print("=" * 60)
    
    update_pvd_files(simulation_dir)
    
    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)
