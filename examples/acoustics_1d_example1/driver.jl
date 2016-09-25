#  Generic driver routine for 1D Clawpack 5.0, using claw1ez
#  Allocation has been moved int claw2ez; this file essentially
#  does nothing, but is being retained because it gives the
#  codebase more flexibility.

include("claw1ez.jl")

function driver()

    claw1ez()

end # program driver

driver()
