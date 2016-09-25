
module SetProb

export
    CParam, read_prob_data

type CParam
    rho   ::Float64
    bulk  ::Float64
    cc    ::Float64
    zz    ::Float64

    beta  ::Float64

    # default constructor
    function CParam()
        new(0.0, 0.0, 0.0, 0.0, 0.0)
    end
end


# Set the material parameters for the acoustic equations

using ReadDataFile

function read_prob_data(parms::CParam)

    r = Reader("setprob.data")

    # density:
    parms.rho  = readflt(r)

    # bulk modulus:
    parms.bulk = readflt(r)

    # sound speed:
    parms.cc = sqrt(parms.bulk/parms.rho)

    # impedance:
    parms.zz = parms.cc*parms.rho

    # beta for initial conditions:
    parms.beta = readflt(r)

    close_reader(r)
end # read_prob_data

end # SetProb

using .SetProb
using Base.Test.@test

function test_read_data()
    const parms = CParam()
    read_prob_data(parms)
    println("CParam = $parms")
    return true
end

# @test test_read_data()


