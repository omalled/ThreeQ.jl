import JLD
import ThreeQ
import PyPlot

@everywhere module Samples

const longembedding = Array{Int64,1}[[0],[5],[3],[6],[2],[4],[12],[9],[13],[10],[14],[11],[15],[23],[16],[21],[18],[22],[17],[20],[28],[24],[30],[26],[29],[27],[31],[39],[35],[36],[34],[38],[32],[37],[45],[40],[44],[42],[47],[55],[49],[54],[51],[53],[48],[52],[60],[58],[62],[56],[61],[59],[63],[71],[65],[70],[67],[69],[64],[68],[76],[72],[77],[75],[78],[73],[79],[87],[82],[86],[80],[85],[81],[84],[92],[89],[95],[88],[94],[90],[93],[91],[187],[188],[184],[190],[186],[191],[185],[189],[181],[179],[182],[176],[183],[177],[180],[172],[168],[174],[170],[175],[169],[173],[165],[162],[166],[161],[167],[163],[164],[156],[153],[158],[152],[159],[155],[157],[149],[144],[151],[146],[150],[145],[148],[140],[136],[142],[138],[143],[137],[141],[133],[128],[135],[129],[134],[130],[132],[124],[122],[126],[123],[127],[121],[125],[117],[113],[118],[114],[119],[112],[116],[108],[105],[111],[104],[110],[106],[109],[101],[97],[102],[98],[100],[96],[192],[197],[195],[198],[193],[199],[194],[196],[204],[201],[205],[202],[206],[214],[211],[213],[208],[215],[209],[212],[220],[216],[222],[217],[221],[218],[223],[231],[225],[230],[224],[229],[226],[228],[236],[234],[238],[233],[237],[232],[239],[247],[242],[245],[241],[246],[240],[244],[252],[251],[253],[250],[255],[249],[254],[262],[256],[260],[257],[261],[269],[266],[270],[264],[268],[265],[271],[279],[275],[278],[274],[277],[272],[276],[284],[280],[285],[283],[286],[282],[378],[380],[376],[381],[373],[371],[375],[369],[374],[368],[372],[364],[360],[366],[363],[367],[361],[365],[357],[353],[358],[352],[359],[355],[356],[348],[346],[350],[347],[351],[344],[349],[341],[336],[342],[338],[343],[339],[340],[332],[329],[335],[328],[333],[325],[322],[326],[320],[327],[321],[324],[316],[312],[318],[313],[319],[314],[317],[309],[307],[311],[305],[310],[304],[308],[300],[296],[302],[298],[303],[299],[301],[293],[291],[295],[290],[294],[288],[384],[389],[385],[391],[387],[390],[386],[388],[396],[395],[398],[394],[397],[393],[399],[407],[400],[405],[401],[406],[403],[404],[412],[410],[414],[408],[413],[409],[415],[423],[416],[421],[419],[422],[418],[420],[428],[425],[430],[426],[429],[424],[431],[439],[432],[437],[433],[438],[434],[436],[444],[440],[445],[442],[446],[441],[447],[455],[449],[454],[451],[453],[448],[452],[460],[458],[463],[471],[465],[469],[467],[470],[464],[468],[476],[472],[478],[474],[477],[473],[479],[475],[571],[572],[569],[575],[568],[574],[570],[573],[565],[560],[567],[563],[566],[561],[564],[556],[554],[559],[553],[558],[555],[557],[549],[545],[551],[546],[550],[544],[548],[540],[537],[543],[538],[542],[536],[541],[533],[528],[534],[530],[535],[531],[532],[524],[520],[526],[522],[527],[521],[525],[517],[514],[518],[512],[519],[513],[516],[508],[507],[511],[505],[510],[506],[509],[501],[499],[503],[497],[502],[496],[500],[492],[488],[494],[490],[495],[489],[493],[485],[483],[486],[482],[484],[480],[576],[581],[579],[583],[578],[582],[577],[580],[588],[584],[589],[586],[590],[587],[591],[599],[592],[598],[594],[597],[595],[596],[604],[600],[606],[602],[605],[603],[607],[615],[609],[614],[608],[613],[621],[618],[622],[617],[623],[631],[627],[629],[625],[630],[626],[628],[636],[632],[637],[634],[638],[633],[639],[647],[643],[646],[642],[645],[641],[644],[652],[648],[654],[650],[653],[651],[655],[663],[659],[661],[658],[662],[656],[660],[668],[665],[761],[767],[760],[766],[763],[765],[757],[753],[759],[754],[758],[752],[756],[748],[744],[750],[746],[751],[745],[749],[741],[736],[742],[738],[743],[737],[740],[732],[728],[734],[729],[735],[731],[733],[725],[723],[726],[720],[724],[716],[713],[719],[712],[717],[714],[718],[710],[705],[708],[700],[698],[703],[699],[701],[693],[691],[695],[690],[694],[688],[692],[684],[680],[686],[681],[687],[683],[685],[677],[674],[678],[675],[679],[673],[676],[672],[768],[773],[771],[775],[769],[774],[770],[772],[780],[778],[782],[779],[781],[776],[783],[791],[786],[790],[787],[789],[785],[788],[796],[794],[798],[792],[797],[795],[799],[807],[800],[806],[803],[805],[813],[811],[814],[810],[812],[808],[815],[823],[831],[839],[833],[838],[835],[837],[832],[836],[844],[842],[846],[841],[845],[840],[847],[855],[849],[854],[850],[853],[848],[852],[860],[859],[862],[856],[863],[857],[861],[858],[954],[956],[952],[958],[953],[957],[949],[944],[951],[945],[950],[946],[948],[940],[936],[942],[937],[943],[939],[941],[933],[930],[934],[929],[935],[931],[932],[924],[920],[926],[921],[927],[923],[925],[917],[915],[918],[914],[919],[913],[916],[908],[904],[910],[905],[911],[907],[909],[901],[896],[902],[897],[903],[898],[900],[892],[891],[895],[889],[894],[890],[893],[885],[883],[887],[882],[886],[880],[884],[876],[872],[878],[875],[879],[873],[877],[869],[866],[870],[865],[871],[867],[868],[864],[960],[967],[961],[966],[963],[965],[962],[964],[972],[971],[974],[969],[975],[983],[979],[982],[978],[981],[977],[980],[988],[984],[990],[985],[989],[986],[991],[999],[992],[997],[993],[998],[994],[996],[1004],[1000],[1005],[1003],[1006],[1001],[1007],[1015],[1010],[1014],[1008],[1013],[1009],[1012],[1020],[1016],[1021],[1018],[1022],[1017],[1023],[1031],[1027],[1029],[1026],[1030],[1025],[1028],[1036],[1035],[1037],[1034],[1038],[1033],[1039],[1047],[1040],[1045],[1042],[1044],[1052],[1050],[1054],[1048],[1055],[1049],[1053],[1051],[1147],[1148],[1145],[1150],[1146],[1151],[1144],[1149],[1141],[1139],[1142],[1138],[1143],[1136],[1140],[1132],[1128],[1135],[1131],[1134],[1129],[1133],[1125],[1120],[1127],[1121],[1126],[1122],[1124],[1116],[1114],[1119],[1113],[1118],[1115],[1117],[1109],[1105],[1111],[1107],[1110],[1106],[1108],[1100],[1096],[1102],[1099],[1103],[1097],[1101],[1093],[1089],[1095],[1091],[1094],[1090],[1092],[1084],[1080],[1086],[1082],[1087],[1081],[1085],[1077],[1074],[1079],[1072],[1078],[1075],[1076],[1068],[1067],[1070],[1066],[1071],[1065],[1069],[1061],[1059],[1060],[1058],[1062],[1057],[1063],[1056]]
klow = 1.0
khigh = 2.0
sigmaobs = 1e-2
#truekb = 0


import ThreeQ
import DataStructures
import PyPlot
#solve 1d groundwater equation assuming the left and right boundary have 0 head, with a spacing of 1 between nodes
function solver(ks::Array{Float64, 1}, rs::Array{Float64, 1})
	I = Array{Int}(3 * (length(ks) - 1) + 2)
	J = Array{Int}(3 * (length(ks) - 1) + 2)
	V = Array{Float64}(3 * (length(ks) - 1) + 2)
	b = Array(Float64, length(ks) + 1)
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
	A = sparse(I, J, V)
	return A \ b
end

function b2f(kb, klow, khigh)
	local ks = Array(Float64, length(kb))
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
	kbs = Array(Int, length(h) - 1, num_reads)
	function finishsolve_helper(m, embans, emb)
		answer = 0.5 * (ThreeQ.DWQMI.unembedanswer(embans["solutions"], emb)' + 1)#convert from ising to qubo
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
	bigkbs = Array(Int, length(h) - 1, num_reads * numsteps)
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
	dwsolver = ThreeQ.DWQMI.getdw2xsys4(Main.mytoken)
	num_reads = 10000
	numcorrect = 0
	srand(0)
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
for a in logspace(0.1, 0.5, 5), b in logspace(-1, -6, 6)
	@time probs[(a, b)] = Samples.getprobabilityoftruekbs(a, b, 100)
	JLD.save("probs.jld", "probs", probs)
end
