SLS_ARM_RATIO = 0.05 # E_1/E0 from (fertella et al. 2025)
GHOST_NODE_SPACING = 8*10**(-6) # 8 micrometers

print("GHOST_NODE_SPACING =", GHOST_NODE_SPACING)
print("SLS_ARM_RATIO =", SLS_ARM_RATIO)
shear_mod = float(input("Enter the shear modulus (Pa): "))
# spacing = float(input("Enter the ghost node spacing (m): "))

E0 = ((shear_mod * 2.9 * GHOST_NODE_SPACING) / SLS_ARM_RATIO)
E1 = SLS_ARM_RATIO * E0

print("E0 =", E0, "Pa")
print("E1 =", E1, "Pa")

# kF = eta_phys / T0
print("-- Sim units --")
E0_sim = 