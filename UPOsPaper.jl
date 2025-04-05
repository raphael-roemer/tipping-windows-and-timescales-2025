using Pkg
Pkg.activate("RoemerAshwin2025")

using DynamicalSystems
using Plots
using LinearAlgebra
using Statistics
using DataFrames
using CSV
using LaTeXStrings
using DelimitedFiles


##############
# find UPOs  #
##############

ds = PredefinedDynamicalSystems.lorenz([0.05,0.05,25.05]; σ = 10.0, ρ = 28.0, β = 8/3)

function periodic_orbit_candidates(ds, total_time;
    sampling_time=0.005, #0.005
    recurrence_threshhold = 0.9, # 0.4
    diagbound = 30, #20
    maxperiod = 1200, # 800
    diagbox_size = 28, #20
    max_line_gap = 50,
    diag_min_xlength = 30, #20
    periodthresh = 15, #9
    semiHausdorffthresh = 2.5 # 3
    )

    T, t = trajectory(ds, total_time; Ttr = 0, Δt = sampling_time)
    Rxyz = RecurrenceMatrix(T, recurrence_threshhold)
    Rixyz = Int8.(Rxyz)
    evsize = length(Rxyz[1,:])-1

    # calculate candidates for periodic orbits
    # We go through the columns of the lower right triangle of the recurrence matrix and always when we find a
    # close recurrence point, we collect all the points from the diagonal "recurrence stripe" in one array. If 
    # the set of points collected in this array has at least some diagonal length, the array will be pushed to 
    # POcand.

    POcand = []
    for i in diagbound:evsize                       # i is the x-axis
        for j in max(1,i-maxperiod):i-diagbound+1   # j is the y-axis
            if Rixyz[j,i] == 1                      # because of matrix-indexing: j (the y-component) comes first
                # if we found one point, we go through the corresponding diagonal
                OneOrbitCandids=[]
                push!(OneOrbitCandids,[i,j])
                k=max(i-5,1)
                line_cont=true
                last_i_added=i
                while k<evsize && line_cont==true
                    for l in max(1,j+(k-i)-diagbox_size):min(k-diagbound+1,j+(k-i)+diagbox_size)
                        if Rixyz[l,k]==1
                            push!(OneOrbitCandids,[k,l])
                            Rixyz[l,k]=0
                            last_i_added = k
                        end
                    end
                    if last_i_added < k-max_line_gap
                        line_cont=false
                    end
                    k+=1
                end
                if last_i_added > diag_min_xlength
                    push!(POcand,OneOrbitCandids)
                end
            end
        end
    end
        
    # now, in each "recurrence stripe" find the representative with the closest recurrence
    POcandids=[]
    for i in 1:length(POcand)
        push!(POcandids,POcand[i][1])
        nor = norm(T[POcand[i][1][1]]-T[POcand[i][1][2]])
        for j in 2:length(POcand[i])
            if norm(T[POcand[i][j][1]]-T[POcand[i][j][2]])<nor
                POcandids[i]=POcand[i][j]
                nor = norm(T[POcand[i][j][1]]-T[POcand[i][j][2]])
            end
        end
    end

    #POperiodl=Array{Int64,1}(undef,Int(2*length(POcandids)))
    POperiodl = []
    for i in 1:length(POcandids)
        push!(POperiodl, Int(POcandids[i][1]-POcandids[i][2])); push!(POperiodl, Int(POcandids[i][1]-POcandids[i][2]))
    end

    # use the starting points POcandids to get the entire trajectory of the POcandids
    # and use the SYMMETRY of the Lorenz System to get all the mirrored orbits (x->-x, y->-y)

    POcandidtrajs = []
    for i in 1:length(POcandids)
        push!(POcandidtrajs, T[POcandids[i][2]:POcandids[i][1]]); 
        push!(POcandidtrajs, T[POcandids[i][2]:POcandids[i][1]])
        for j in 1:length(POcandidtrajs[end])
            POcandidtrajs[end][j]=SVector{3}(-POcandidtrajs[end][j][1],-POcandidtrajs[end][j][2],POcandidtrajs[end][j][3])
        end
        #println(i," ", i+1, " ", POcandidtrajs[i+1]==POcandidtrajs[i])
    end

    # Now, identify all POs that are very similar AND have a period that is n∈N times the period of each other
    # To achieve this, follow these steps:
    # - take all pairs of POs with similar period or roughly n∈N times the period of each other. 
    # - For each pair, go through all the points of both the two POs and check whether the semi-Hausdorff distance 
    #   of each of these points to the other PO exceeds a threshold. 
    # - If it exceeds the threshhold --> the POs are not identified
    #   If it does not exceed the threshhold --> the POs are identified and the one with the smallest period is stored

    identifiedPOs=[]
    for i in 1:length(POcandidtrajs)
        for j in i+1:length(POcandidtrajs)
            smallper=min(POperiodl[i],POperiodl[j])
            largeper=max(POperiodl[i],POperiodl[j])
            # check whether the longer period is roughly n∈N times the period of the shorter
            if min(abs(mod(largeper,smallper)-0),abs(mod(largeper,smallper)-smallper))<periodthresh
                # if yes, go through many of the points of both the two POs and check whether the semi-Hausdorff distance 
                #   of each of these points to the other PO exceeds a threshold. 
                distij = zeros(length(1:2:POperiodl[i]),length(1:2:POperiodl[j]))
                for k in 1:length(1:2:POperiodl[i])
                    for l in 1:length(1:2:POperiodl[j])
                        distij[k,l] = norm(POcandidtrajs[i][(1:2:POperiodl[i])[k]]-POcandidtrajs[j][(1:2:POperiodl[j])[l]])
                        #distij[k,l] = norm(T[(POcandids[i][2]:2:POcandids[i][1])[k]]-T[(POcandids[j][2]:2:POcandids[j][1])[l]])
                    end
                end
                
                kmins=zeros(length(distij[:,1]))
                for k in 1:length(distij[:,1])
                    kmins[k] = minimum(distij[k,:])
                end
                maxdist_i_to_j = maximum(kmins)
                #println("Max dist of a point in ",i," to a point in ",j," is: ",round(maxdist_i_to_j,digits=2))
                
                lmins=zeros(length(distij[1,:]))
                for l in 1:length(distij[1,:])
                    lmins[l] = minimum(distij[:,l])
                end
                #println(round.(lmins,digits=4))
                maxdist_j_to_i = maximum(lmins)
                #println("Max dist of a point in ",j," to a point in ",i," is: ",round(maxdist_j_to_i,digits=2))


                if max(maxdist_i_to_j,maxdist_j_to_i) < semiHausdorffthresh
                    if POperiodl[i] > POperiodl[j]
                        push!(identifiedPOs,i)
                    else
                        push!(identifiedPOs,j)

                    end
                end
            end
        end
    end

    POcandidates=[]
    POs_length=[]
    for i in 1:length(POcandidtrajs)
        if i ∉ identifiedPOs
            push!(POcandidates,POcandidtrajs[i])
            push!(POs_length,POperiodl[i])
        end
    end

    A = zeros(length(POcandidates),2)
    A[:,1]=POs_length
    A[:,2]=1:length(POcandidates)
    POslength_sortl_temp=sortslices(A, dims=1)
    POs_sortl=[]
    for i in 1:length(POcandidates)
        push!(POs_sortl,POcandidates[Int(POslength_sortl_temp[i,2])])
    end
    POslength_sortl=Int.(POs_length[Int.(POslength_sortl_temp[:,2])])

    # sort the POs such that of every rotated pair the PO with the larger xmean comes first.
    po_stats = periodic_orbits_stats(POs_sortl)
    POs_sortl2=POs_sortl[1:end]
    POslength_sortl2=POslength_sortl[1:end]
    for i in 1:length(POs_sortl)-1
        if POslength_sortl[i] == POslength_sortl[i+1] && po_stats[i,1] < po_stats[i+1,1]
            println(i)
            POs_sortl2[i] = POs_sortl[i+1]
            POs_sortl2[i+1] = POs_sortl[i]
            POslength_sortl2[i] = POslength_sortl[i+1]
            POslength_sortl2[i+1] = POslength_sortl[i]
        end
    end

    return POs_sortl2, POslength_sortl2, T, t, Rixyz, identifiedPOs, POcandidtrajs
end
function periodic_orbits_stats(POss)
    po_stats=zeros(length(POss),4)
    for i in 1:length(POss)
        xmean = mean(POss[i][:,1])
        recurrence = norm(POss[i][1,:]-POss[i][length(POss[i]),:])
        xmax = maximum(POss[i][:,1])
        xmin = minimum(POss[i][:,1])
        po_stats[i,:]=[xmean,recurrence,xmax,xmin]
    end
    po_stats
end

POst, POslt, Traj, Trajt, Rixyz, identifiedPOs, POcandidtrajs = periodic_orbit_candidates(ds,100); # 300
POst
POslt
POs=POst[1:20]
POsl=POslt[1:20]
POs[16:20]=POst[18:22]
POsl[16:20]=POslt[18:22]

# print a recurrence matrix to illustrate it in the paper
T, t = trajectory(ds, 20; Ttr = 0, Δt = 0.005)
Rxyz = RecurrenceMatrix(T, 0.9)
Rintxyz = Int8.(Rxyz)
Mrecpl=heatmap(Rintxyz[1:4000,1:4000], size=(500,500),  color=cgrad(:devon, 2, categorical = true, rev=true),cbar=false, xlabel=L"i", ylabel=L"j")


# plot all POs
PO_p=[];
col_count=1;
for candcounter in 1:length(POs)
    line_style_c=:solid
    markerstyle_c=:none
    col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
    if candcounter > 1 && length(POs[candcounter])!=length(POs[candcounter-1])
        col_count+=1
        col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
        if candcounter < length(POs)-1 && length(POs[candcounter])!=length(POs[candcounter+1]) # symmetric POs are dashed, mirrored ones are both solid line
            line_style_c=:dash
        end
    end
    if candcounter > 1 && length(POs[candcounter])==length(POs[candcounter-1]) # symmetric POs are dashed, mirrored ones are both solid line
        line_style_c=:dot
        markerstyle_c=:circle
        #println(candcounter)
    end
    push!(PO_p,plot(POs[candcounter][:,1],POs[candcounter][:,2],POs[candcounter][:,3],
        size=(750,600),
        xlabel=L"x",
        ylabel=L"y",
        zlabel=L"z",
        label=false,  #L"PO \ %$candcounter",
        title=L"PO \ %$candcounter",
        xlim=(-18,18),
        ylim=(-24,24),
        zlim=(5,45),
        linewidth=2,
        color=col_c,
        linestyle= line_style_c,
        markershape=markerstyle_c,
        markerstrokecolor=col_c,
        markersize=2,
        xguidefontsize=20,yguidefontsize=20,zguidefontsize=20,
        xtickfontsize=14,ytickfontsize=14,ztickfontsize=14,
        titlefont=(20,"Computer Modern")
        ))
end;
popl=plot(
    PO_p[1], 
    PO_p[2], 
    PO_p[3], 
    PO_p[4],
    PO_p[5], 
    PO_p[6], 
    PO_p[7], 
    PO_p[8],
    PO_p[9], 
    PO_p[10], 
    PO_p[11], 
    PO_p[12],
    PO_p[13], 
    PO_p[14], 
    PO_p[15], 
    PO_p[16],
    PO_p[17],
    PO_p[18],
    PO_p[19],
    PO_p[20],
    layout=(5,4), size=(2300,3000)
)
# plot single PO
#candcounter=32
#plot(POst[candcounter][:,1],POst[candcounter][:,2],POst[candcounter][:,3], size=(1000,800), linestyle=:dot)


po_stats = periodic_orbits_stats(POs)
# histogram(po_stats[:,1], xlabel="xmean", bins=30)
# histogram(po_stats[:,2], xlabel="recurrence", bins=30)
# histogram(po_stats[:,3], xlabel="xmax", bins=30)
# histogram(po_stats[:,4], xlabel="xmin", bins=30)
# histogram(POsl, xlabel="duration", bins=30)
# plot(POsl,log.(1:length(POsl)))

function POs_y1_to_long_trajectories(POs, POsl, Traj, leng;
    timestep_long=0.05
    )

    POs_y_long=[]
    for i in 1:length(POs)
        PO_y_long=zeros(round(Int64,leng/timestep_long)+1)
        for j in 1:round(Int64,leng/timestep_long)+1
            #200 =1/0.005  (indices per timestep)
            PO_y_long[j]=POs[i][Int64(floor(mod((1/0.005)*((j-1)*timestep_long),length(POs[i]))))+1,1]
        end
        push!(POs_y_long,PO_y_long)
    end
    return POs_y_long, 0:timestep_long:leng
end
POs_y1_long, POs_y1_l_long = POs_y1_to_long_trajectories(POs, POsl, Traj,10000, timestep_long=0.05)