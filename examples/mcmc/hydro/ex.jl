import JLD
import ThreeQ
import PyPlot

@everywhere module Samples

const longembedding = Array{Int64,1}[[0],[5],[3],[6],[2],[4],[12],[9],[14],[10],[15],[11],[13],[21],[19],[22],[17],[20],[16],[23],[31],[24],[29],[27],[30],[25],[28],[36],[34],[38],[33],[39],[32],[37],[45],[40],[44],[42],[47],[55],[49],[52],[48],[54],[51],[53],[61],[56],[62],[59],[63],[57],[60],[68],[64],[69],[67],[70],[65],[71],[79],[73],[78],[75],[76],[72],[77],[85],[81],[86],[82],[87],[80],[84],[92],[89],[95],[88],[94],[90],[93],[91],[187],[189],[184],[190],[186],[191],[185],[188],[180],[177],[183],[176],[182],[179],[181],[173],[169],[175],[170],[174],[168],[172],[164],[163],[165],[162],[167],[161],[166],[158],[155],[156],[153],[159],[152],[157],[149],[144],[151],[146],[150],[145],[148],[140],[136],[142],[138],[143],[137],[141],[133],[130],[135],[129],[132],[131],[134],[126],[123],[125],[121],[127],[122],[124],[116],[112],[119],[114],[118],[113],[117],[109],[106],[110],[104],[111],[105],[108],[100],[97],[102],[99],[101],[96],[192],[197],[195],[198],[193],[199],[194],[196],[204],[202],[206],[203],[205],[213],[210],[214],[208],[212],[211],[215],[223],[218],[220],[216],[222],[217],[221],[229],[227],[231],[224],[230],[225],[228],[236],[234],[238],[233],[237],[232],[239],[247],[242],[244],[241],[245],[240],[246],[254],[248],[253],[250],[255],[249],[252],[260],[257],[261],[256],[262],[270],[266],[271],[265],[268],[264],[269],[277],[275],[279],[274],[278],[272],[276],[284],[280],[285],[283],[286],[282],[378],[381],[377],[380],[372],[369],[374],[371],[375],[368],[373],[365],[360],[366],[363],[367],[361],[364],[356],[355],[357],[352],[359],[353],[358],[350],[344],[351],[346],[348],[347],[349],[341],[336],[342],[338],[343],[339],[340],[332],[329],[335],[328],[333],[325],[323],[324],[321],[327],[320],[326],[318],[313],[319],[314],[317],[312],[316],[308],[304],[309],[306],[311],[305],[310],[302],[296],[300],[299],[303],[297],[301],[293],[291],[295],[290],[294],[288],[384],[389],[385],[391],[387],[390],[386],[388],[396],[392],[398],[394],[399],[393],[397],[405],[403],[406],[401],[404],[400],[407],[415],[409],[413],[408],[414],[410],[412],[420],[416],[422],[418],[423],[419],[421],[429],[426],[430],[425],[428],[424],[431],[439],[432],[436],[435],[438],[434],[437],[445],[441],[447],[440],[446],[442],[444],[452],[450],[453],[451],[454],[449],[455],[463],[458],[460],[468],[466],[471],[464],[469],[465],[470],[478],[474],[476],[472],[477],[473],[479],[475],[571],[573],[570],[574],[568],[575],[569],[572],[564],[561],[566],[563],[567],[560],[565],[557],[555],[558],[553],[559],[554],[556],[548],[544],[551],[545],[549],[547],[550],[542],[536],[540],[537],[543],[538],[541],[533],[528],[534],[530],[535],[531],[532],[524],[520],[525],[523],[527],[521],[526],[518],[513],[516],[515],[519],[512],[517],[509],[506],[510],[505],[511],[507],[508],[500],[496],[502],[497],[503],[499],[501],[493],[489],[495],[490],[494],[488],[492],[484],[482],[486],[481],[485],[480],[576],[580],[577],[583],[579],[582],[578],[581],[589],[586],[590],[587],[591],[584],[588],[596],[595],[597],[594],[598],[592],[599],[607],[600],[606],[614],[608],[615],[609],[612],[604],[603],[605],[613],[621],[618],[622],[617],[623],[631],[627],[628],[624],[630],[625],[629],[637],[634],[638],[633],[639],[635],[636],[644],[641],[645],[642],[646],[643],[647],[655],[651],[653],[650],[654],[648],[652],[660],[659],[662],[656],[663],[658],[661],[669],[665],[761],[767],[760],[765],[763],[764],[756],[752],[758],[754],[759],[753],[757],[749],[745],[751],[746],[750],[744],[748],[740],[737],[743],[738],[741],[739],[742],[734],[730],[732],[729],[735],[728],[733],[725],[723],[726],[720],[724],[716],[713],[719],[712],[717],[715],[718],[710],[705],[709],[701],[699],[703],[698],[700],[692],[688],[694],[690],[695],[691],[693],[685],[683],[687],[681],[686],[680],[684],[676],[674],[678],[675],[679],[673],[677],[672],[768],[772],[769],[774],[770],[775],[771],[773],[781],[778],[782],[779],[783],[776],[780],[788],[785],[789],[787],[790],[786],[791],[799],[795],[797],[792],[796],[793],[798],[806],[800],[807],[803],[805],[813],[811],[814],[810],[812],[808],[815],[823],[831],[825],[830],[826],[829],[827],[828],[836],[832],[837],[835],[838],[833],[839],[847],[840],[845],[841],[846],[842],[844],[852],[848],[854],[850],[855],[849],[853],[861],[859],[860],[856],[863],[857],[862],[858],[954],[957],[953],[958],[952],[956],[948],[944],[951],[945],[950],[946],[949],[941],[939],[943],[937],[942],[936],[940],[932],[931],[935],[929],[933],[930],[934],[926],[921],[927],[920],[924],[922],[925],[917],[915],[918],[914],[919],[913],[916],[908],[904],[911],[907],[909],[906],[910],[902],[897],[903],[898],[900],[899],[901],[893],[890],[894],[889],[895],[891],[892],[884],[880],[886],[882],[887],[883],[885],[877],[873],[879],[875],[878],[872],[876],[868],[866],[870],[865],[871],[867],[869],[864],[960],[967],[961],[966],[963],[965],[962],[964],[972],[968],[975],[969],[974],[982],[977],[981],[978],[980],[979],[983],[991],[986],[989],[985],[990],[984],[988],[996],[994],[999],[992],[998],[993],[997],[1005],[1000],[1004],[1003],[1006],[1001],[1007],[1015],[1010],[1014],[1008],[1012],[1009],[1013],[1021],[1017],[1023],[1018],[1022],[1016],[1020],[1028],[1025],[1030],[1026],[1029],[1027],[1031],[1039],[1032],[1037],[1034],[1038],[1033],[1036],[1044],[1040],[1045],[1042],[1047],[1055],[1048],[1054],[1050],[1052],[1049],[1053],[1051],[1147],[1151],[1146],[1149],[1144],[1150],[1145],[1148],[1140],[1136],[1143],[1138],[1142],[1139],[1141],[1133],[1129],[1134],[1131],[1135],[1128],[1132],[1124],[1122],[1127],[1120],[1125],[1121],[1126],[1118],[1112],[1119],[1114],[1116],[1113],[1117],[1109],[1105],[1111],[1107],[1110],[1106],[1108],[1100],[1096],[1102],[1099],[1103],[1097],[1101],[1093],[1089],[1095],[1088],[1092],[1090],[1094],[1086],[1081],[1087],[1083],[1085],[1080],[1084],[1076],[1075],[1078],[1072],[1079],[1074],[1077],[1069],[1065],[1071],[1066],[1070],[1067],[1068],[1060],[1058],[1062],[1059],[1063],[1057],[1061],[1056]]
klow = 1.0
khigh = 2.0
sigmaobs = 1e-2
#truekb = 0


import ThreeQ
import DataStructures
import PyPlot
import Random
import SparseArrays
#solve 1d groundwater equation assuming the left and right boundary have 0 head, with a spacing of 1 between nodes
function solver(ks::Array{Float64, 1}, rs::Array{Float64, 1})
	I = Array{Int}(undef, 3 * (length(ks) - 1) + 2)
	J = Array{Int}(undef, 3 * (length(ks) - 1) + 2)
	V = Array{Float64}(undef, 3 * (length(ks) - 1) + 2)
	b = Array{Float64}(undef, length(ks) + 1)
	I[1] = 1; J[1] = 1; V[1] = 1.0; b[1] = length(ks)
	I[2] = length(ks) + 1; J[2] = length(ks) + 1; V[2] = 1.0; b[2] = 0.0
	push!(I, length(ks) + 1); push!(J, length(ks) + 1); push!(V, 1.); b[end] = 0.0
	l = 3
	dx = 1
	for i = 2:length(ks)
		b[i] = -rs[i - 1]
		I[l] = i
		J[l] = i - 1
		V[l] = ks[i - 1] / dx^2
		I[l + 1] = i
		J[l + 1] = i + 1
		V[l + 1] = ks[i] / dx^2
		I[l + 2] = i
		J[l + 2] = i
		V[l + 2] = -(ks[i - 1] + ks[i]) / dx^2
		l += 3
	end
	A = SparseArrays.sparse(I, J, V)
	return A \ b
end

function b2f(kb, klow, khigh)
	local ks = Array{Float64}(undef, length(kb))
	for i = 1:length(kb)
		if kb[i] == 1
			ks[i] = khigh
		else
			ks[i] = klow
		end
	end
	return ks
end

function solver(kb, rs, klow, khigh)
	ks = b2f(kb, klow, khigh)
	return solver(ks, rs)
end

function sample(h, klow, khigh, rs; num_reads=100, solver=ThreeQ.DWQMI.defaultsolver, kwargs...)
	Q = zeros(length(h) - 1, length(h) - 1)
    for i = 2:length(h) - 1
		constterm = klow * h[i + 1] - klow * h[i] - klow * h[i] + klow * h[i - 1] + rs[i - 1]
		qicoeff = (khigh - klow) * (h[i + 1] - h[i])
		qim1coeff = (khigh - klow) * (h[i - 1] - h[i])
		Q[i, i] += qicoeff ^ 2 + 2 * constterm * qicoeff
		Q[i - 1, i - 1] += qim1coeff ^ 2 + 2 * constterm * qim1coeff
		Q[i, i - 1] += 2 * qicoeff * qim1coeff
    end
	#=
	fig, axs = PyPlot.subplots(2, 3, figsize=(16, 9))
	axs[1][:plot](sort(hi), "k.")
	axs[2][:plot](sort(collect(values(J))), "k.")
	axs[3][:plot](h)
	axs[3][:plot](Samples.solver(fill(Float64(klow), length(q.value)), rs))
	=#
	embeddings = longembedding[1:size(Q, 1)]
	adjacency = ThreeQ.DWQMI.getadjacency(solver)
	kbs = Array{Int}(undef, length(h) - 1, num_reads)
	function finishsolve_helper(m, embans, emb)
		answer = 0.5 * (ThreeQ.DWQMI.unembedanswer(embans["solutions"], emb)' .+ 1)#convert from ising to qubo
		i = 1
		for j = 1:size(answer, 2)
			occurrences = embans["num_occurrences"][j]
			#=
			if j <= 1
				@show map(Int, answer[:, j]) == truekb
				@show occurrences
				axs[4][:plot](Samples.solver(klow + (khigh - klow) * answer[:, j], rs))
				axs[5][:plot](h - Samples.solver(klow + (khigh - klow) * answer[:, j], rs))
				axs[6][:plot](h - Samples.solver(klow + (khigh - klow) * rand([0.0, 1.0], size(answer, 1)), rs))
			end
			=#
			for l = 1:occurrences
				kbs[:, i] = answer[:, j]
				i += 1
			end
		end
	end
    ThreeQ.solvesapi!(Q; auto_scale=true, solver=solver, num_reads=num_reads, embeddings=embeddings, finishsolve_helper=finishsolve_helper, kwargs...)
	#=
	display(fig); println()
	PyPlot.close(fig)
	=#
	return kbs
end

function sample(h, klow, khigh, rs, myloglikelihood; burnin=100, numsteps=100, num_reads=100, kwargs...)
	kbs = sample(h, klow, khigh, rs; num_reads=num_reads, kwargs...)
	bigkbs = Array{Int}(undef, length(h) - 1, num_reads * numsteps)
	i = 1
	for j = 1:size(kbs, 2)
		bigkbs[:, i:i + numsteps - 1], _ = ThreeQ.sample(kbs[:, j], burnin, numsteps, myloglikelihood)
		i += numsteps
    end
    return bigkbs
end

function makellhood(numhs)
	global truekb = rand([0, 1], numhs - 1)
	truers = zeros(numhs - 2)
	trueh = solver(truekb, truers, klow, khigh)
	obsh = trueh + sigmaobs * randn(length(trueh))
	@show sigmaobs
	@show klow, khigh
	function myloglikelihood(kb)
		h = solver(kb, truers, klow, khigh)
		return -sum((h - obsh).^2) / (2 * sigmaobs^2)
	end
	return myloglikelihood, truekb, truers, obsh
end

function exactcompare(numhs, numcomparisons, burnin, numsteps)
	myloglikelihood, truekb, truers, obsh = makellhood(numhs)
	dwsolver = ThreeQ.DWQMI.getdw2xsys4(Main.mytoken)
	trueprobs = Dict(zip(ThreeQ.enumerateprobabilities(length(truekb), myloglikelihood)...))
	classicals = SharedArray{Float64}(numcomparisons)
	dwaves = SharedArray{Float64}(numcomparisons)
	for j = 1:numcomparisons
		#do the classical sampling
		classicalsamples, _ = ThreeQ.sample(rand([0, 1], length(truekb)), burnin, numsteps, myloglikelihood)
		classicalmcmcprobs = DataStructures.DefaultDict{Array{eltype(classicalsamples), 1}, Float64}(0.0)
		for i = 1:size(classicalsamples, 2)
			classicalmcmcprobs[classicalsamples[:, i]] += 1 / size(classicalsamples, 2)
		end
		classicals[j] = ThreeQ.bhatdist(trueprobs, classicalmcmcprobs)
	end
	#do the d-wave sampling
	alldwavesamples = sample(obsh, klow, khigh, truers, myloglikelihood; num_reads=numcomparisons, solver=dwsolver, burnin=burnin, numsteps=numsteps, token=Main.mytoken, timeout=Inf)
	for j = 1:numcomparisons
		dwavesamples = alldwavesamples[:, 1 + (j - 1) * numsteps:j * numsteps]
		dwavemcmcprobs = DataStructures.DefaultDict{Array{eltype(dwavesamples), 1}, Float64}(0.0)
		for i = 1:size(dwavesamples, 2)
			dwavemcmcprobs[dwavesamples[:, i]] += 1 / size(dwavesamples, 2)
		end
		dwaves[j] = ThreeQ.bhatdist(trueprobs, dwavemcmcprobs)
	end
	return dwaves, classicals
end

function countsamples(numhs, maxsteps, numreads)
	myloglikelihood, truekb, truers, obsh = makellhood(numhs)
	dwsolver = ThreeQ.DWQMI.getdw2xsys4(Main.mytoken)
	kbs = sample(obsh, klow, khigh, truers; num_reads=numreads, solver=dwsolver, token=Main.mytoken, timeout=Inf)
	targetllhood = myloglikelihood(kbs[:, 1])
	samples, llhoods = ThreeQ.sample(rand([0, 1], length(truekb)), 0, maxsteps, myloglikelihood)
	bestsofar = copy(llhoods)
	for i = 2:length(llhoods)
		bestsofar[i] = max(bestsofar[i - 1], llhoods[i])
	end
	return targetllhood, bestsofar
end

function getprobabilityoftruekbs(thiskhigh, thissigmaobs, numrepeats)
	global khigh = thiskhigh
	global sigmaobs = thissigmaobs
	#dwsolver = ThreeQ.DWQMI.getdw2xsys4(Main.mytoken)
	dwsolver = ThreeQ.DWQMI.getdw2xvfyc(Main.mytoken)
	num_reads = 100
	numcorrect = 0
	Random.seed!(0)
	for j = 1:numrepeats
		myloglikelihood, truekb, truers, obsh = makellhood(length(longembedding) + 1)
		kbs = sample(obsh, klow, khigh, truers; num_reads=num_reads, solver=dwsolver, token=Main.mytoken, timeout=Inf)
		for i = 1:size(kbs, 2)
			if kbs[:, i] == truekb
				numcorrect += 1
			end
		end
	end
	return numcorrect / (num_reads * numrepeats)
end

end#end Samples Module

#dw, cl = Samples.exactcompare(20, 100, 10^2, 10^2)
#tllh, bestsofar = Samples.countsamples(900, 10^4, 10000)
probs = Dict()
#for a in logspace(0.1, 0.5, 5), b in logspace(-1, -6, 6)
for a in 10.0 .^ range(0.1, 0.5; length=2), b in 10.0 .^ range(-1, -6; length=2)
	@show a, b
	@time probs[(a, b)] = Samples.getprobabilityoftruekbs(a, b, 10)
	JLD.save("probs_new.jld", "probs", probs)
end
