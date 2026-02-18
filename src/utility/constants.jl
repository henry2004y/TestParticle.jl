using PhysicalConstants.CODATA2018
using Unitful: ustrip, @u_str

const qₑ = -ustrip(u"C", ElementaryCharge)  # electron charge, [C]
const mₑ = ustrip(u"kg", ElectronMass)      # electron mass, [kg]
const qᵢ = ustrip(u"C", ElementaryCharge)   # proton charge, [C]
const mᵢ = ustrip(u"kg", ProtonMass)        # proton mass, [kg]
const c = ustrip(u"m/s", SpeedOfLightInVacuum)       # speed of light, [m/s]
const μ₀ = ustrip(u"H/m", VacuumMagneticPermeability)          # Vacuum permeability, [H/m]
const ϵ₀ = ustrip(u"F/m", VacuumElectricPermittivity)       # Vacuum permittivity, [F/m]
const kB = ustrip(u"J/K", BoltzmannConstant)   # Boltzmann constant, [m²kg/(s²K)]
const Rₑ = 6371.0e3           # Earth radius, [m]
const BMoment_Earth = SVector(0.0, 0.0, -7.94e22) # [V*s/(A*m)]

const Proton = Species(mᵢ, qᵢ)
const Electron = Species(mₑ, qₑ)

Ion(m, q = 1) = Species(m * mᵢ, q * qᵢ)
Ion(; m = 1, q = 1) = Species(m * mᵢ, q * qᵢ)
