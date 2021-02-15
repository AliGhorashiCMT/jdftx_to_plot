filename = (@__DIR__)*"/../data/reference_txt.csv"
open(filename; write=true, create=true, truncate=true) do io
	write(io, "a reference")
end

