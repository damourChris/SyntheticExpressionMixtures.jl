# Script to setup the conda environment 
using CondaPkg
using Preferences
using Libdl
using UUIDs

# Setup the channels
CondaPkg.resolve()

target_rhome = joinpath(CondaPkg.envdir(), "lib", "R")
if Sys.iswindows()
    target_libr = joinpath(target_rhome, "bin", Sys.WORD_SIZE == 64 ? "x64" : "i386",
                           "R.dll")
else
    target_libr = joinpath(target_rhome, "lib", "libR.$(Libdl.dlext)")
end

const RCALL_UUID = UUID("6f49c342-dc21-5d91-9882-a32aef131414")

# Check if LocalPreferences.toml exists
if !isfile(joinpath(dirname(@__DIR__), "LocalPreferences.toml"))
    # Create a LocalPreferences.toml file
    open(joinpath(dirname(@__DIR__), "LocalPreferences.toml"), "w") do io
        return write(io, """
               [RCall]
               Rhome = "$target_rhome"
               libR = "$target_libr"
               """)
    end
else
    set_preferences!(RCALL_UUID, "Rhome" => target_rhome, "libR" => target_libr;
                     force=true)
end