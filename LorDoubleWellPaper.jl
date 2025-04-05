

# Roemer, Ashwin 2025
#########################################################################################
#########################################################################################
###  Make sure to run the file "UPOs.jl" first to compute the unstable periodic orbits
###  needed for the following code!
#########################################################################################
#########################################################################################



###############################################################################
# Plot of Lorenz Double Well ensembles, non-tipping, partial tipping, tipping #
###############################################################################
# run the system for three different pairs of (a, η) use different Lorenz initial 
# conditions for each pair of (a, η), but always the same initial condition for the Double Well
mutable struct DWParameters_fullSyst
    γ::Float64
    μ::Float64
    a::Float64
end
function DWfullSystem(x, p, t)
    σ = 10.0; ρ = 28.0; β = 8/3;

    dx1 = σ*(x[2] - x[1])
    dx2 = x[1]*(ρ - x[3]) - x[2]
    dx3 = x[1]*x[2] -β*x[3]
    dx4 = (-x[4]^3+3x[4] + p.μ + p.a*x[1])/p.γ

    return SVector{4}(dx1, dx2, dx3, dx4)
end
function DWtypical_trajs_plotvals(μs, a, tmax, ensemble_size;
    γ=1,
    Δt=0.05,
    )

    ret_traj=zeros(ensemble_size*length(μs), length(0.0:Δt:tmax))

    inits = [0.1,0.1,25.1,-1.5]
    for i in 1:length(μs)
        p=DWParameters_fullSyst(γ,μs[i],a)
        ds = ContinuousDynamicalSystem(DWfullSystem, inits, p, t0=0)
        step!(ds,tmax,true)
        inits=[get_state(ds)[1],get_state(ds)[2],get_state(ds)[3],-1.5]
        for j in 1:ensemble_size
            ds2 = ContinuousDynamicalSystem(DWfullSystem, inits, p, t0=0)

            ret_traj[(i-1)*ensemble_size+j,:] = trajectory(ds2,tmax,Δt=Δt)[1][:,4]
            inits=[get_state(ds2)[1],get_state(ds2)[2],get_state(ds2)[3],-1.5]
        end
    end

    pl=plot((0.0:Δt:tmax)./γ,[ret_traj[k,:] for k ∈ 1:ensemble_size],label=false, color="blue", xlabel=L"t", ylabel=L"x")
    plot!(pl,(0.0:Δt:tmax)./γ,[ret_traj[k,:] for k ∈ ensemble_size+1:2*ensemble_size],label=false, color="orange")
    plot!(pl,(0.0:Δt:tmax)./γ,[ret_traj[k,:] for k ∈ 2*ensemble_size+1:3*ensemble_size],label=false, color="lime")
end
DWtypical_trajs_plotvals([1.7,1.85,2.4], 0.03, 150, 7)


################
# a vs η plots #
################
# Lorenz-Double-Well system forced by UPOs
mutable struct DWParameters_PO
    γ::Float64
    μ::Float64
    a::Float64
    POs_y1_long::Vector
    timestep_UPOs::Float64
end
function DW_UPOs(x, p, t)
    dx1 = (-x[1]^3+3x[1] + p.μ + p.a*p.POs_y1_long[round(Int,t/p.timestep_UPOs)+1])/p.γ

    return SVector{1}(dx1)
end
# Lorenz-Double-Well system forced by constant value (e.g. the mean or maximum x-value of the PO)
mutable struct DWParameters_const
    γ::Float64
    μ::Float64
    a::Float64
    forcing_val::Float64
end
function DW_const(x, p, t)
    dx1 = (-x[1]^3+3x[1] + p.μ + p.a*p.forcing_val)/p.γ

    return SVector{1}(dx1)
end


function thresh_cross_a_mu(PO_nrs,fixed_forcing_vect;
    γ=1,
    as=0.0:0.2:0.4,
    μ_1_bounds=(0.0,2.0),
    Δt=2,
    accuracy=5e-3
    )

    mu_1_a=[]
    for i in 1:length(PO_nrs)
        println("PO ", PO_nrs[i])
        mu_1_a_1PO = zeros(length(as),2)
        for j in 1:length(as)
            (μ_n1, μ_n2) = μ_1_bounds
            while abs(μ_n2-μ_n1) > accuracy
                μ_temp = (μ_n2+μ_n1)/2
                p=DWParameters_const(γ, μ_temp, as[j],fixed_forcing_vect[i])
                ds = ContinuousDynamicalSystem(DW_const, [-1.5], p, t0=0)
                
                tipped=false
                while tipped==false && current_time(ds) < 150
                    step!(ds,Δt,true)
                    if get_state(ds)[1] > 0.0 
                        tipped = true
                        (μ_n1, μ_n2) = (μ_n1,μ_temp)
                    end
                end
                if tipped==false # not tipped
                    (μ_n1, μ_n2) = (μ_temp,μ_n2)
                end
            end
            mu_1_a_1PO[j,:] = [(μ_n2+μ_n1)/2, as[j]]
        end
        push!(mu_1_a,mu_1_a_1PO)
    end
    return mu_1_a
end
function thresh_cross_a_mu(PO_nrs;
    γ=1,
    as=0.0:0.2:0.4,
    μ_1_bounds=(0.0,2.0),
    Δt=0.5, #0.05,
    accuracy=5e-3,
    timestepUPOs=0.05  # should be <= to Δt
    )


    POs_y1_long, POs_y1_l_long = POs_y1_to_long_trajectories(POs, POsl, Traj,10000, timestep_long=timestepUPOs)

    mu_1_a=[]
    for i in 1:length(PO_nrs)
        # println("PO ", PO_nrs[i])
        mu_1_a_1PO = zeros(length(as),2)
        for j in 1:length(as)
            println("PO ", PO_nrs[i], ", a = ",as[j])
            (μ_n1, μ_n2) = μ_1_bounds
            while abs(μ_n2-μ_n1) > accuracy
                μ_temp = (μ_n2+μ_n1)/2
                p=DWParameters_PO(γ, μ_temp, as[j],POs_y1_long[i],timestepUPOs)
                ds = ContinuousDynamicalSystem(DW_UPOs, [-1.5], p, t0=0)

                tipped=false
                while tipped==false && current_time(ds) < 150
                    step!(ds,Δt,true)
                    if get_state(ds)[1] > 0.0 
                        tipped = true
                        (μ_n1, μ_n2) = (μ_n1,μ_temp)
                    end
                end
                if tipped==false # not tipped
                    (μ_n1, μ_n2) = (μ_temp,μ_n2)
                end
            end
            mu_1_a_1PO[j,:] = [(μ_n2+μ_n1)/2, as[j]]
        end
        push!(mu_1_a,mu_1_a_1PO)
    end
    return mu_1_a
end
function makePlot_mu_a(mu_a_data,gamma)
    P1=plot(xlabel=L"\eta", ylabel=L"a",title=L"\gamma = %$(gamma)");
    col_count=1
    for pocounter in 1:length(POs)
        line_style_c=:solid
        markerstyle_c=:none
        col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
        if pocounter > 1 && length(POs[pocounter])!=length(POs[pocounter-1])
            col_count+=1
            col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
            if pocounter < length(POs)-1 && length(POs[pocounter])!=length(POs[pocounter+1]) # symmetric POs are dashed, mirrored ones are both solid line
                line_style_c=:dash
            end
        end
        if pocounter > 1 && length(POs[pocounter])==length(POs[pocounter-1]) # symmetric POs are dashed, mirrored ones are both solid line
            line_style_c=:dot
            markerstyle_c=:circle
        end
        plot!(P1,mu_a_data[pocounter][:,1],mu_a_data[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            xlim=(1.,2.5),
            linewidth=2,
            color=col_c,
            linestyle= line_style_c,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=2,
            xguidefontsize=20,yguidefontsize=20,
            xtickfontsize=12,ytickfontsize=12,
            titlefont=(20,"Computer Modern"),
            legendfontsize=12
            )
    end
    return(P1)
end

POs_cands=1:length(POs)
as=0.0:0.01:0.04   # as=0.0:0.005:0.04
μ_1_bounds=(1.2,2.5)
# γ = 0 :
# The Lorenz-POs run much slower than the DW system, and thus, we expect the maximum 
# of the POs' forcing to be relevant
mu_1_a_0_gamma=thresh_cross_a_mu(POs_cands, po_stats[:,3]; γ=1, as=as, μ_1_bounds=μ_1_bounds)
makePlot_mu_a(mu_1_a_0_gamma,0)
# γ = 0.01:
mu_1_a_gamma0_01=thresh_cross_a_mu(POs_cands,γ=0.01,Δt=0.05, as=as, μ_1_bounds=μ_1_bounds,timestepUPOs=0.05)
#makePlot_mu_a(mu_1_a_gamma0_01,0.01)
# γ = 0.1 :
mu_1_a_gamma0_1=thresh_cross_a_mu(POs_cands,γ=0.1,Δt=0.05, as=as, μ_1_bounds=μ_1_bounds,timestepUPOs=0.05)
#makePlot_mu_a(mu_1_a_gamma0_1,0.1)
# γ = 0.5:
mu_1_a_gamma0_5=thresh_cross_a_mu(POs_cands,γ=0.5,Δt=0.05, as=as, μ_1_bounds=μ_1_bounds,timestepUPOs=0.05)
#makePlot_mu_a(mu_1_a_gamma0_5,0.5)
# γ = 1 :
mu_1_a=thresh_cross_a_mu(POs_cands,γ=1, Δt=0.05, as=as, μ_1_bounds=μ_1_bounds,timestepUPOs=0.05)
#makePlot_mu_a(mu_1_a,1)
# γ = 10 :
mu_1_a_gamma10=thresh_cross_a_mu(POs_cands,γ=10,Δt=0.1, as=as, μ_1_bounds=μ_1_bounds,timestepUPOs=0.1)
#makePlot_mu_a(mu_1_a_gamma10,10)
# infinite γ
# The Lorenz-POs run much faster than the DW system, and thus, the Stommel system only 
# "sees" the mean of the POs' forcing.
mu_1_a_inf_gamma=thresh_cross_a_mu(POs_cands, po_stats[:,1]; γ=1, as=as, μ_1_bounds=μ_1_bounds)
makePlot_mu_a(mu_1_a_inf_gamma,"∞")





######################################
# plot of μ vs γ for fixed a=0.08  #
######################################
function makeplot_mu_gamma_fixed_a()
    plotvals = zeros(length(POs),6,2)
    for i in 1:length(POs)
        plotvals[i,1,:] = [mu_1_a_0_gamma[i][length(as),1],0.008]  #[eta_vals,γ_vals]
        plotvals[i,2,:] = [mu_1_a_gamma0_01[i][length(as),1],0.01]
        plotvals[i,3,:] = [mu_1_a_gamma0_1[i][length(as),1],0.1]
        plotvals[i,4,:] = [mu_1_a[i][length(as),1],1.0]
        plotvals[i,5,:] = [mu_1_a_gamma10[i][length(as),1],10.0]
        plotvals[i,6,:] = [mu_1_a_inf_gamma[i][length(as),1],12.0]
    end
    mu_a_at_0_08=plot(xlabel=L"\eta", ylabel=L"\gamma",title=L"a = 0.04");
    col_count=1
    for pocounter in 1:length(POs)
        line_style_c=:solid
        markerstyle_c=:none
        col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
        if pocounter > 1 && length(POs[pocounter])!=length(POs[pocounter-1])
            col_count+=1
            col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
            if pocounter < length(POs)-1 && length(POs[pocounter])!=length(POs[pocounter+1]) # symmetric POs are dashed, mirrored ones are both solid line
                line_style_c=:dash
            end
        end
        if pocounter > 1 && length(POs[pocounter])==length(POs[pocounter-1]) # symmetric POs are dashed, mirrored ones are both solid line
            line_style_c=:dot
            markerstyle_c=:circle
        end
        plot!(plotvals[pocounter,2:end-1,1],plotvals[pocounter,2:end-1,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            yscale=:log10,
            #xlim=(1.05,1.35),
            #ylim=(-40,40),
            #yticks=plotvals[pocounter,1:end-1,2],
            linewidth=2,
            color=col_c,
            linestyle= line_style_c,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=2,
            xguidefontsize=20,yguidefontsize=20,
            xtickfontsize=12,ytickfontsize=12,
            titlefont=(20,"Computer Modern"),
            legendfontsize=12
            )

        plot!([plotvals[pocounter,[1,6],1]],[plotvals[pocounter,[1,6],2]],
            size=(750,600),
            label=false,
            alpha=1,
            #yscale=:ln,
            color=col_c,
            seriestype=:scatter,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=4)
    end
    return mu_a_at_0_08
end
mu_a_at_0_08pl=makeplot_mu_gamma_fixed_a()







#############################################################
# γ vs. a of the the intersection of two PO1 and another PO #
#############################################################
# we use the relation between a and η given by PO1
μPO1(a)= 2. - a*po_stats[1,3]

# now, go through several γs and do a bisection in a with η given by ηPO1(a) 
function DWthresh_cross_a_gamma(PO_nr;
    aLowerUpper=(0,0.23),
    γs=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8],
    Δt=0.05, #0.05,
    accuracy=5e-3,
    finaltime=150,
    timestepUPOs=0.05  # should be <= to Δt
    )

    POs_y1_long, POs_y1_l_long = POs_y1_to_long_trajectories(POs, POsl, Traj,10000, timestep_long=timestepUPOs)
    
    a_gamma=zeros(length(γs),2)
    for (γ,γ_) in enumerate(γs)
        println(L"\gamma ", γ_)
        (a_n1, a_n2) = aLowerUpper
        while abs(a_n2-a_n1) > accuracy
            a_temp = (a_n2+a_n1)/2
            p=DWParameters_PO(γ_, μPO1(a_temp), a_temp,POs_y1_long[PO_nr],timestepUPOs)
            ds = ContinuousDynamicalSystem(DW_UPOs, [-1.5], p, t0=0)
            
            tipped=false
            while tipped==false && current_time(ds) < finaltime
                step!(ds,Δt,true)
                if get_state(ds)[1] > 0.0 
                    tipped = true
                    (a_n1, a_n2) = (a_n1,a_temp)
                end
            end
            if tipped==false # not tipped
                (a_n1, a_n2) = (a_temp,a_n2)
            end
        end
        a_gamma[γ,:] = [γ_, (a_n1+a_n2)/2]
    end
    return a_gamma
end
a_gamma=DWthresh_cross_a_gamma(4,γs=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1])
plot(a_gamma[1:end,1],a_gamma[1:end,2],xlabel=L"\gamma",ylabel=L"a",label=false,xguidefontsize=20,yguidefontsize=20,xtickfontsize=12,ytickfontsize=12,titlefont=(20,"Computer Modern"),guidefontsize=20);
plot!(a_gamma[1:end,1],(a_gamma[1:end,1].^2).*21,xlabel=L"\gamma",ylabel=L"a",label=false)

plot(a_gamma[1:end,1],a_gamma[1:end,2],xlabel=L"\gamma",ylabel=L"a",yscale=:log10,xscale=:log10,label=false,xguidefontsize=20,yguidefontsize=20,xtickfontsize=12,ytickfontsize=12,titlefont=(20,"Computer Modern"));
plot!(a_gamma[1:end,1],(a_gamma[1:end,1].^2).*21,xlabel=L"\gamma",ylabel=L"a",yscale=:log10,xscale=:log10,label=false)





###################################
# section 3) Plot non-UPO forcing #
###################################
# for each eta and a, use a single new orbit starting at the final point of the previous orbit. 
# Fix the total integration time and then do the plot again with a longer integration time


# Lorenz-Double well system forced by general Lorenz trajectory
# store time when it tipped
function yellow_blue_plotvals_time(μs, as, tmax;
    γ=1,
    Δt=0.05,
    t_interm=5
    )

    inits = [0.1,0.1,25.1,-1.5]
    tipped_aμ=zeros(length(μs),length(as))
    for i in 1:length(μs)
        for j in 1:length(as)

            # run  Lorenz for t_interm Lorenz-timesteps to decorelater subsequent runs
            p = DWParameters_fullSyst(γ,μs[i],as[j])
            dst = ContinuousDynamicalSystem(DWfullSystem, inits, p, t0=0)
            traj,trajt = trajectory(dst, t_interm; Ttr = t_interm-0.01)

            inits = [traj[end][1],traj[end][2],traj[end][3],-1.5]
            ds = ContinuousDynamicalSystem(DWfullSystem, inits, p, t0=0)
            while tipped_aμ[i,j]==0.0 && current_time(ds) < tmax
                step!(ds,Δt,true)
                if get_state(ds)[4] > 0.0 
                    tipped_aμ[i,j] = current_time(ds)
                end
            end
            inits = [get_state(ds)[1],get_state(ds)[2],get_state(ds)[3],-1.5]
        end
    end
    return tipped_aμ
end
muvals = 1.25:0.01:3.0
avals = 0.0:0.0005:0.04
yb_plotvals_g0_01_tm150=yellow_blue_plotvals_time(muvals, avals, 100, γ=0.01,Δt=0.001)
yb_plotvals_g0_1_tm150=yellow_blue_plotvals_time(muvals, avals, 150, γ=0.1,Δt=0.005)
yb_plotvals_g0_5_tm150=yellow_blue_plotvals_time(muvals, avals, 150, γ=0.5,Δt=0.05)
yb_plotvals_g1_tm150=yellow_blue_plotvals_time(muvals, avals, 150, γ=1)
yb_plotvals_g10_tm150=yellow_blue_plotvals_time(muvals, avals, 150, γ=10,Δt=0.1)

# sort now into four time intervals:
function sortybplotvals4(yb_vals,upper,mid,lower)
    outp = zeros(length(yb_vals[:,1]),length(yb_vals[1,:]))
    for i in 1:length(yb_vals[:,1])
        for j in 1:length(yb_vals[1,:])
            if upper > yb_vals[i,j] > mid
                outp[i,j]=1.
            elseif mid > yb_vals[i,j] > lower
                outp[i,j]=2.
            elseif lower > yb_vals[i,j]>0
                outp[i,j]=3.
            end
        end
    end
    outp
end
cat_plotvals0_01=sortybplotvals4(yb_plotvals_g0_01_tm150,150,0.05,0.0175)
cat_plotvals0_1=sortybplotvals4(yb_plotvals_g0_1_tm150,150,0.5,0.12)
cat_plotvals0_5=sortybplotvals4(yb_plotvals_g0_5_tm150,150,1.0,0.5)
cat_plotvals1=sortybplotvals4(yb_plotvals_g1_tm150,150,3,1.25)
cat_plotvals10=sortybplotvals4(yb_plotvals_g10_tm150,150,20,12)

####################################
# combine the heatmap with PO plots
####################################
# extend heatmap values to the size needed
mustep=round((muvals[end]-muvals[1])/length(muvals),digits=3)
ext_muvals=0.4:mustep:3.05
function extend_yb_plotvals(yb_plotvals,mustep,ext_muvals)
    ext_yb_plotvals=zeros(length(avals),length(ext_muvals))
    ext_yb_plotvals[:,1:length(0.4:mustep:muvals[1])].=0.0
    ext_yb_plotvals[:,length(0.4:mustep:muvals[1])+1 : length(0.4:mustep:muvals[1])+length(muvals)].=yb_plotvals'
    ext_yb_plotvals[:,length(0.4:mustep:muvals[1])+1+length(muvals):end].=3.0
    
    return ext_yb_plotvals
end

function make_overlayed_plot_overlayed0_01(hm_vals,mustep,ext_muvals)
    ext_hm_vals = extend_yb_plotvals(hm_vals,mustep,ext_muvals)

    P1_gamma0_01=plot(xlabel=L"\eta", ylabel=L"a",title=L"\gamma = 0.01");
    heatmap!(P1_gamma0_01,ext_muvals, avals, ext_hm_vals,colorbar=false,color=cgrad(:binary,rev=true),xguidefontsize=20,yguidefontsize=20,xtickfontsize=12,ytickfontsize=12,titlefont=(20,"Computer Modern"));
    col_count=1
    for pocounter in 1:length(POs)
        line_style_c=:solid
        markerstyle_c=:none
        col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
        if pocounter > 1 && length(POs[pocounter])!=length(POs[pocounter-1])
            col_count+=1
            col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
            if pocounter < length(POs)-1 && length(POs[pocounter])!=length(POs[pocounter+1]) # symmetric POs are dashed, mirrored ones are both solid line
                line_style_c=:dash
            end
        end
        if pocounter > 1 && length(POs[pocounter])==length(POs[pocounter-1]) # symmetric POs are dashed, mirrored ones are both solid line
            line_style_c=:dot
            markerstyle_c=:circle
        end
        plot!(P1_gamma0_01,mu_1_a_gamma0_01[pocounter][:,1],mu_1_a_gamma0_01[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            #xlim=(0.4,1.85),# (1.028,1.083),  
            #ylim=(0.0,0.081),# (0.0325,0.046),
            linewidth=2,
            color=col_c,
            linestyle= line_style_c,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=2,
            xguidefontsize=20,yguidefontsize=20,
            xtickfontsize=12,ytickfontsize=12,
            titlefont=(20,"Computer Modern"),
            legendfontsize=12)
    end
    return P1_gamma0_01
end
make_overlayed_plot_overlayed0_01(cat_plotvals0_01,mustep,ext_muvals)

function make_overlayed_plot_overlayed0_1(hm_vals,mustep,ext_muvals)
    ext_hm_vals = extend_yb_plotvals(hm_vals,mustep,ext_muvals)

    P1_gamma0_1=plot(xlabel=L"\eta", ylabel=L"a",title=L"\gamma = 0.1");
    heatmap!(P1_gamma0_1,ext_muvals, avals, ext_hm_vals,colorbar=false,color=cgrad(:binary,rev=true),alpha=1.);
    col_count=1
    for pocounter in 1:length(POs)
        line_style_c=:solid
        markerstyle_c=:none
        col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
        if pocounter > 1 && length(POs[pocounter])!=length(POs[pocounter-1])
            col_count+=1
            col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
            if pocounter < length(POs)-1 && length(POs[pocounter])!=length(POs[pocounter+1]) # symmetric POs are dashed, mirrored ones are both solid line
                line_style_c=:dash
            end
        end
        if pocounter > 1 && length(POs[pocounter])==length(POs[pocounter-1]) # symmetric POs are dashed, mirrored ones are both solid line
            line_style_c=:dot
            markerstyle_c=:circle
        end
        plot!(P1_gamma0_1,mu_1_a_gamma0_1[pocounter][:,1],mu_1_a_gamma0_1[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            #xlim=(0.4,1.85),# (1.028,1.083),  
            #ylim=(0.0,0.081),# (0.0325,0.046),
            linewidth=2,
            color=col_c,
            linestyle= line_style_c,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=2,
            xguidefontsize=20,yguidefontsize=20,
            xtickfontsize=12,ytickfontsize=12,
            titlefont=(20,"Computer Modern"),
            legendfontsize=12)
    end
    return P1_gamma0_1
end
make_overlayed_plot_overlayed0_1(cat_plotvals0_1,mustep,ext_muvals)

function make_overlayed_plot_overlayed0_5(hm_vals,mustep,ext_muvals)
    ext_hm_vals = extend_yb_plotvals(hm_vals,mustep,ext_muvals)

    P1_gamma0_5=plot(xlabel=L"\eta", ylabel=L"a",title=L"\gamma = 0.5");
    heatmap!(P1_gamma0_5,ext_muvals, avals, ext_hm_vals,colorbar=false,color=cgrad(:binary,rev=true),alpha=1.);
    col_count=1
    for pocounter in 1:length(POs)
        line_style_c=:solid
        markerstyle_c=:none
        col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
        if pocounter > 1 && length(POs[pocounter])!=length(POs[pocounter-1])
            col_count+=1
            col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
            if pocounter < length(POs)-1 && length(POs[pocounter])!=length(POs[pocounter+1]) # symmetric POs are dashed, mirrored ones are both solid line
                line_style_c=:dash
            end
        end
        if pocounter > 1 && length(POs[pocounter])==length(POs[pocounter-1]) # symmetric POs are dashed, mirrored ones are both solid line
            line_style_c=:dot
            markerstyle_c=:circle
        end
        plot!(P1_gamma0_5,mu_1_a_gamma0_5[pocounter][:,1],mu_1_a_gamma0_5[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            #xlim=(0.4,1.85),# (1.028,1.083),  
            #ylim=(0.0,0.081),# (0.0325,0.046),
            linewidth=2,
            color=col_c,
            linestyle= line_style_c,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=2,
            xguidefontsize=20,yguidefontsize=20,
            xtickfontsize=12,ytickfontsize=12,
            titlefont=(20,"Computer Modern"),
            legendfontsize=12)
    end
    return P1_gamma0_5
end
make_overlayed_plot_overlayed0_5(cat_plotvals0_5,mustep,ext_muvals)

function make_overlayed_plot_overlayed1(hm_vals,mustep,ext_muvals)
    ext_hm_vals = extend_yb_plotvals(hm_vals,mustep,ext_muvals)

    P1_gamma1=plot(xlabel=L"\eta", ylabel=L"a",title=L"\gamma = 1");
    heatmap!(P1_gamma1,ext_muvals, avals, ext_hm_vals,colorbar=false,color=cgrad(:binary,rev=true),alpha=1.);
    col_count=1
    for pocounter in 1:length(POs)
        line_style_c=:solid
        markerstyle_c=:none
        col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
        if pocounter > 1 && length(POs[pocounter])!=length(POs[pocounter-1])
            col_count+=1
            col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
            if pocounter < length(POs)-1 && length(POs[pocounter])!=length(POs[pocounter+1]) # symmetric POs are dashed, mirrored ones are both solid line
                line_style_c=:dash
            end
        end
        if pocounter > 1 && length(POs[pocounter])==length(POs[pocounter-1]) # symmetric POs are dashed, mirrored ones are both solid line
            line_style_c=:dot
            markerstyle_c=:circle
        end
        plot!(P1_gamma1,mu_1_a[pocounter][:,1],mu_1_a[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            #xlim=(0.4,1.85),# (1.028,1.083),  
            #ylim=(0.0,0.081),# (0.0325,0.046),
            linewidth=2,
            color=col_c,
            linestyle= line_style_c,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=2,
            xguidefontsize=20,yguidefontsize=20,
            xtickfontsize=12,ytickfontsize=12,
            titlefont=(20,"Computer Modern"),
            legendfontsize=12)
    end
    return P1_gamma1
end
make_overlayed_plot_overlayed1(cat_plotvals1,mustep,ext_muvals)

function make_overlayed_plot_overlayed10(hm_vals,mustep,ext_muvals)
    ext_hm_vals = extend_yb_plotvals(hm_vals,mustep,ext_muvals)

    P1_gamma10=plot(xlabel=L"\eta", ylabel=L"a",title=L"\gamma = 10");
    heatmap!(P1_gamma10,ext_muvals, avals, ext_hm_vals,colorbar=false,color=cgrad(:binary,rev=true),alpha=1.);
    col_count=1
    for pocounter in 1:length(POs)
        line_style_c=:solid
        markerstyle_c=:none
        col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
        if pocounter > 1 && length(POs[pocounter])!=length(POs[pocounter-1])
            col_count+=1
            col_c=cgrad(:gist_rainbow, 11, categorical = true)[col_count]
            if pocounter < length(POs)-1 && length(POs[pocounter])!=length(POs[pocounter+1]) # symmetric POs are dashed, mirrored ones are both solid line
                line_style_c=:dash
            end
        end
        if pocounter > 1 && length(POs[pocounter])==length(POs[pocounter-1]) # symmetric POs are dashed, mirrored ones are both solid line
            line_style_c=:dot
            markerstyle_c=:circle
        end
        plot!(P1_gamma10,mu_1_a_gamma10[pocounter][:,1],mu_1_a_gamma10[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            #xlim=(0.4,1.85),# (1.028,1.083),  
            #ylim=(0.0,0.081),# (0.0325,0.046),
            linewidth=2,
            color=col_c,
            linestyle= line_style_c,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=2,
            xguidefontsize=20,yguidefontsize=20,
            xtickfontsize=12,ytickfontsize=12,
            titlefont=(20,"Computer Modern"),
            legendfontsize=12)
    end
    return P1_gamma10
end
make_overlayed_plot_overlayed10(cat_plotvals10,mustep,ext_muvals)