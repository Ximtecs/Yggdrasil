using Printf
using LinearAlgebra
using Base.Threads

include("SnapshotMetaStructs.jl")
include("PrintMetaStructs.jl")
include("ValueParsing.jl")
include("SnapshotParser.jl")
include("IDXConversing.jl")
include("HelperFunctions.jl")
include("FindPatches.jl")
include("CartesianGrid.jl")
include("Scaling.jl")

include("LoadPatch.jl")
include("LoadSnapshot.jl")
include("LoadSnapshotParticles.jl")
include("LoadMultipleSnapshots.jl")
include("LoadMultipleSnapshotsParticles.jl")

