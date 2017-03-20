import Requests

function uint82bits!(a, x, n)
	a[n + 8] = (x & 0x01) != 0x00
	a[n + 7] = (x & 0x02) != 0x00
	a[n + 6] = (x & 0x04) != 0x00
	a[n + 5] = (x & 0x08) != 0x00
	a[n + 4] = (x & 0x10) != 0x00
	a[n + 3] = (x & 0x20) != 0x00
	a[n + 2] = (x & 0x40) != 0x00
	a[n + 1] = (x & 0x80) != 0x00
end

function numbytespersolution(numbits)
	x = div(numbits, 8)
	y = mod(numbits, 8)
	return y == 0 ? x : x + 1
end

function getanswer(problem, token, url)
	result = Requests.json(Requests.get("$url/problems/$problem"; headers=Dict("X-Auth-Token"=>token), tls_conf=Requests.TLS_NOVERIFY))
	if !haskey(result, "answer")
		error("error in the JSON result:\n$result")
	end
	answer = result["answer"]
	active_variables = reinterpret(Int32, base64decode(answer["active_variables"]))
	energies = reinterpret(Float64, base64decode(answer["energies"]))
	num_occurrences = reinterpret(Int32, base64decode(answer["num_occurrences"]))
	uint8solutions = base64decode(answer["solutions"])
	solutions1d = Array(Bool, 8 * length(uint8solutions))
	n = 0
	for i = 1:length(uint8solutions)
		uint82bits!(solutions1d, uint8solutions[i], n)
		n += 8
	end
	numbytes = numbytespersolution(length(active_variables))
	numsolutions = div(length(uint8solutions), numbytes)
	solutions = fill(3, numsolutions, answer["num_variables"])
	if result["type"] == "qubo"
		for i = 1:numsolutions, j = 1:length(active_variables)
			solutions[i, active_variables[j] + 1] = solutions1d[(i - 1) * numbytes * 8 + j] ? 1 : 0
		end
	elseif result["type"] == "ising"
		for i = 1:numsolutions, j = 1:length(active_variables)
			solutions[i, active_variables[j] + 1] = solutions1d[(i - 1) * numbytes * 8 + j] ? 1 : -1
		end
	else
		error("unknown problem type: $(result["type"])")
	end
	answer["active_variables"] = active_variables
	answer["energies"] = energies
	answer["num_occurrences"] = num_occurrences
	answer["solutions"] = solutions
	return answer
end

function savebqpjson(Q, filename; description="", id=0, metadata=Dict(), variable_ids=collect(0:size(Q, 1)-1), variable_domain="boolean", scale=1.0, offset=0.0)
	linear_terms = []
	quadratic_terms = []
	for i = 1:size(Q, 1)
		if Q[i, i] != 0
			push!(linear_terms, Dict("id"=>variable_ids[i], "coeff"=>Q[i, i]))
		end
		for j = 1:i - 1
			if Q[i, j] + Q[j, i] != 0
				push!(quadratic_terms, Dict("id_tail"=>variable_ids[j], "id_head"=>variable_ids[i], "coeff"=>Q[i, j] + Q[j, i]))
			end
		end
	end
	bqpdict = Dict("version"=>1.0.0, "id"=>id, "metadata"=>metadata, "variable_ids"=>variable_ids, "variable_domain"=>variable_domain, "scale"=>scale, "offset"=>offset, "linear_terms"=>linear_terms, "quadratic_terms"=>quadratic_terms)
	f = open(filename, "w")
	write(f, JSON.json(bqpdict))
	close(f)
end
