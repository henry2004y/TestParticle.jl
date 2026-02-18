using PhysicalConstants.CODATA2018
using Unitful: ustrip

const qₑ = -ustrip(ElementaryCharge) # electron charge, [C]
const mₑ = ustrip(ElectronMass) # electron mass, [kg]
const qᵢ = ustrip(ElementaryCharge) # proton charge, [C]
const mᵢ = ustrip(ProtonMass) # proton mass, [kg]
const c = ustrip(SpeedOfLightInVacuum) # speed of light, [m/s]
const μ₀ = ustrip(VacuumMagneticPermeability) # Vacuum permeability, [H/m]
const ϵ₀ = ustrip(VacuumElectricPermittivity) # Vacuum permittivity, [F/m]
const kB = ustrip(BoltzmannConstant) # Boltzmann constant, [J/K]
const Rₑ = 6371.0e3 # Earth radius, [m]
const BMoment_Earth = SVector(0.0, 0.0, -7.94e22) # [V*s/(A*m)]

const Proton = Species(mᵢ, qᵢ)
const Electron = Species(mₑ, qₑ)

Ion(m, q = 1) = Species(m * mᵢ, q * qᵢ)
Ion(; m = 1, q = 1) = Species(m * mᵢ, q * qᵢ)
