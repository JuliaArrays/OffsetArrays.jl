module ReadDataFile

export Reader, readstr, readbool, readint, readflt, readintarray, readfltarray, close_reader

type Reader
       
    fname::String
    istream::IOStream
    line::String
    lineno::Int

    # a constructor: opens a file, sets lineno
    function Reader(fname::String)
        istream = open(fname, "r")
        r = new(fname,istream, "", 0)
        finalizer(r, close_reader)
        return r
    end
end

function readln(r::Reader)
    const skip_rexp = r"^\s*(?:#|$)"
    const rexp = r"^\s*(.*?)\s*(?:(#|=:)\s*(.*?)\s*$|$)"

    while true
        r.line = readline(r.istream)
        r.lineno += 1
        # skip empty lines and comments starting from #
        if match(skip_rexp, r.line) == nothing
            break
        end
    end
    m = match(rexp, r.line)
    r.line = m.captures[1]
end

function readstr(r::Reader)
    readln(r)
    return r.line
end

function readbool(r::Reader)
    readln(r)
    if r.line != "T" && r.line != "F"
        local fname = r.fname, lineno = r.lineno
        error("$fname($lineno): error parsing boolean")
    end
    return r.line == "T"
end

function readint(r::Reader)
    readln(r)
    return parse(Int,r.line)
end

function readflt(r::Reader)
    readln(r)
    return parse(Float64, r.line)
end

function readraw(r::Reader)
    readln(r)
    return split(r.line)
end

function readarray(r::Reader, parser, len)
    buf = readraw(r)
    buflen = length(buf)
    @assert(buflen == len)
    return [parser(buf[i]) for i = 1:buflen]
end

readintarray(r::Reader, len::Int) = readarray(r, x -> parse(Int, x), len)
readfltarray(r::Reader, len::Int) = readarray(r, x -> parse(Float64, x), len)

function close_reader(r::Reader)
    Base.close(r.istream)
end

# ```jlcon
# julia> include("ReadDataFile.jl")
# 
# julia> using ReadDataFile
# 
# julia> r = Reader("claw.data")
# -I- claw.data 6 blank lines were skipped
# Reader("claw.data",IOStream(<file claw.data>),"1",6,r"^\s*(.*?)\s*(?:#\s*(.*?)\s*$|$)")
# 
# julia> readint(r)
# 1
# 
# julia> readflt(r)
# -1.0
# ```

end # ReadDataFile
