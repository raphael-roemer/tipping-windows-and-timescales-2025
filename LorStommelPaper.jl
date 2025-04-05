

# Roemer, Ashwin 2025
#########################################################################################
#########################################################################################
###  Make sure to run the file "UPOs.jl" first to compute the unstable periodic orbits
###  needed for the following code!
#########################################################################################
#########################################################################################



###########################################################################
# Plot of lorenz Stommel ensembles, non-tipping, partial tipping, tipping #
###########################################################################
# run the system for three different pairs of (a, η) use different Lorenz initial 
# conditions for each pair of (a, η), but always the same initial condition for the Stommel
mutable struct Parameters_fullSyst
    γ::Float64
    ξ_1::Float64
    η_1::Float64
    ζ::Float64
    a::Float64
end
function fullSystem(x, p, t)
    σ = 10.0; ρ = 28.0; β = 8/3;

    dx1 = σ*(x[2] - x[1])
    dx2 = x[1]*(ρ - x[3]) - x[2]
    dx3 = x[1]*x[2] -β*x[3]
    dx4 = ((p.ξ_1+p.a*x[1])-x[4]*(1+abs(x[4]-x[5])))/p.γ
    dx5 = ((p.η_1+p.a*x[1])-x[5]*(p.ζ+abs(x[4]-x[5])))/p.γ

    return SVector{5}(dx1, dx2, dx3, dx4, dx5)
end
function typical_trajs_plotvals(ηs, a, tmax, ensemble_size;
    γ=1,
    Δt=0.05,
    )

    Ψ_dLorenzStommel_traj=zeros(ensemble_size*length(ηs), length(0.0:Δt:tmax))

    inits = [0.1,0.1,25.1,1.7,1.0]
    for i in 1:length(ηs)
        p=Parameters_fullSyst(γ,3, ηs[i], 0.3, a)
        dLorenzStommel = ContinuousDynamicalSystem(fullSystem, inits, p, t0=0)
        step!(dLorenzStommel,tmax)
        inits=[get_state(dLorenzStommel)[1],get_state(dLorenzStommel)[2],get_state(dLorenzStommel)[3],1.7,1.0]
        for j in 1:ensemble_size
            dLorenzStommel2 = ContinuousDynamicalSystem(fullSystem, inits, p, t0=0)
            dLorenzStommel_traj2=trajectory(dLorenzStommel2,tmax,Δt=0.05)
            for k in 1:length(dLorenzStommel_traj2[1][:])
                Ψ_dLorenzStommel_traj[(i-1)*ensemble_size+j,k]=dLorenzStommel_traj2[1][k][4]-dLorenzStommel_traj2[1][k][5]
            end
            inits=[get_state(dLorenzStommel2)[1],get_state(dLorenzStommel2)[2],get_state(dLorenzStommel2)[3],1.7,1.0]
            #reinit!(dLorenzStommel, [get_state(dLorenzStommel)[1],get_state(dLorenzStommel)[2],get_state(dLorenzStommel)[3],1.7,1.0])
        end
    end

    pl=plot(0.0:0.05:150.0,[Ψ_dLorenzStommel_traj[k,:] for k ∈ 1:ensemble_size],label=false, color="blue", xlabel=L"t", ylabel=L"\Psi")
    plot!(pl,0.0:0.05:150.0,[Ψ_dLorenzStommel_traj[k,:] for k ∈ ensemble_size+1:2*ensemble_size],label=false, color=:orange)
    plot!(pl,0.0:0.05:150.0,[Ψ_dLorenzStommel_traj[k,:] for k ∈ 2*ensemble_size+1:3*ensemble_size],label=false, color="lime")

    return pl
end
LSyb_plotvals=typical_trajs_plotvals([1.0,1.175,1.5], 0.04, 150, 7)





################
# a vs η plots #
################
# Lorenz-Stommel system forced by POs
mutable struct Params1
    γ::Float64
    PO_nr::Int64
    ξ_1::Float64
    η_1::Float64
    ζ::Float64
    a::Float64
    PO_y1::Vector
    timestepUPOs::Float64
end
function G_POs(x, p, t)
    dx1 = ((p.ξ_1+p.a*p.PO_y1[round(Int,t/p.timestepUPOs)+1])-x[1]*(1+abs(x[1]-x[2])))/p.γ
    dx2 = ((p.η_1+p.a*p.PO_y1[round(Int,t/p.timestepUPOs)+1])-x[2]*(p.ζ+abs(x[1]-x[2])))/p.γ

    return SVector{2}(dx1, dx2)
end

# Lorenz-Stommel system forced by constant value (e.g. the mean or maximum x-value of the PO)
mutable struct Parameters_const_forcin1
    γ::Float64
    forcing_vals::Float64
    ξ_1::Float64
    η_1::Float64
    ζ::Float64
    a::Float64
end
function G_const(x, p, t)
    dx1 = ((p.ξ_1+p.a*p.forcing_vals)-x[1]*(1+abs(x[1]-x[2])))/p.γ
    dx2 = ((p.η_1+p.a*p.forcing_vals)-x[2]*(p.ζ+abs(x[1]-x[2])))/p.γ

    return SVector{2}(dx1, dx2) #SVector has a better performance for less then 100 entries
end


function thresh_cross_a_eta(PO_nrs,fixed_forcing_vect;
    γ=1,
    as=0.0:0.2:0.4,
    η_1_bounds=(0.0,2.0),
    Δt=2,
    accuracy=5e-3
    )

    eta_1_a=[]
    for i in 1:length(PO_nrs)
        println("PO ", PO_nrs[i])
        eta_1_a_1PO = zeros(length(as),2)
        for j in 1:length(as)
            (η_n1, η_n2) = η_1_bounds
            while abs(η_n2-η_n1) > accuracy
                η_temp = (η_n2+η_n1)/2
                p_const_forc=Parameters_const_forcin1(γ,fixed_forcing_vect[i],3, η_temp, 0.3, as[j])
                dStommel_const = ContinuousDynamicalSystem(G_const, [1.7,1.0], p_const_forc, t0=0)
                
                tipped=false
                while tipped==false && current_time(dStommel_const) < 150
                    step!(dStommel_const,Δt)
                    Ψ = get_state(dStommel_const)[1]-get_state(dStommel_const)[2]
                    if Ψ < 0.1 
                        tipped = true
                        (η_n1, η_n2) = (η_n1,η_temp)
                    end
                end
                if tipped==false # not tipped
                    (η_n1, η_n2) = (η_temp,η_n2)
                end
            end
            eta_1_a_1PO[j,:] = [(η_n2+η_n1)/2, as[j]]
        end
        push!(eta_1_a,eta_1_a_1PO)
    end
    return eta_1_a
end
function thresh_cross_a_eta(PO_nrs;
    γ=1,
    as=0.0:0.2:0.4,
    η_1_bounds=(0.0,2.0),
    Δt=0.05, #0.05,
    accuracy=5e-3,
    timestepUPOs=0.05,  # should be <= to Δt
    tmax=150
    )

    POs_y1_long, POs_y1_l_long = POs_y1_to_long_trajectories(POs, POsl, Traj,10000, timestep_long=timestepUPOs)

    eta_1_a=[]
    for i in 1:length(PO_nrs)
        println("PO ", PO_nrs[i])
        eta_1_a_1PO = zeros(length(as),2)
        for j in 1:length(as)
            (η_n1, η_n2) = η_1_bounds
            while abs(η_n2-η_n1) > accuracy
                η_temp = (η_n2+η_n1)/2
                p=Params1(γ,PO_nrs[i],3, η_temp, 0.3, as[j],POs_y1_long[i],timestepUPOs)
                dStommel = ContinuousDynamicalSystem(G_POs, [1.7,1.0], p, t0=0)
                
                tipped=false
                while tipped==false && current_time(dStommel) < tmax
                    step!(dStommel,Δt,true)
                    Ψ = get_state(dStommel)[1]-get_state(dStommel)[2]
                    if Ψ < 0.1 
                        tipped = true
                        (η_n1, η_n2) = (η_n1,η_temp)
                    end
                end
                if tipped==false # not tipped
                    (η_n1, η_n2) = (η_temp,η_n2)
                end
            end
            eta_1_a_1PO[j,:] = [(η_n2+η_n1)/2, as[j]]
        end
        push!(eta_1_a,eta_1_a_1PO)
    end
    return eta_1_a
end
function makePlot_eta_a(eta_a_data,gamma)
    P1=plot(xlabel=L"\eta_1", ylabel=L"a",title=L"\gamma = %$(gamma)");
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
        plot!(P1,eta_a_data[pocounter][:,1],eta_a_data[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            xlim=(0.4,1.6),
            #ylim=(-40,40),
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
    return P1
end

POs_cands=1:length(POs)
as_LS=0.0:0.005:0.08
η_1_bounds_LS=(0.0,3.0)

# γ = 0 :
# The Lorenz-POs run much slower than the Stommel system, and thus, we expect the maximum 
# of the POs' forcing to be relevant
eta_1_a_0_gamma=thresh_cross_a_eta(POs_cands, po_stats[:,3]; γ=1, as=as_LS, η_1_bounds=η_1_bounds_LS)
makePlot_eta_a(eta_1_a_0_gamma,0)

# γ = 0.01 :
eta_1_a_gamma0_01=thresh_cross_a_eta(POs_cands,γ=0.01,Δt=0.05, as=as_LS, η_1_bounds=η_1_bounds_LS,timestepUPOs=0.05,tmax=50)
#makePlot_eta_a(eta_1_a_gamma0_01,0.01)
# γ = 0.05 :
eta_1_a_gamma0_05=thresh_cross_a_eta(POs_cands,γ=0.05, as=as_LS, η_1_bounds=η_1_bounds_LS,timestepUPOs=0.05,tmax=50)
#makePlot_eta_a(eta_1_a_gamma0_05,0.05)
# γ = 0.1 :
eta_1_a_gamma0_1=thresh_cross_a_eta(POs_cands,γ=0.1,Δt=0.05, as=as_LS, η_1_bounds=η_1_bounds_LS,timestepUPOs=0.05)
#makePlot_eta_a(eta_1_a_gamma0_1,0.1)
# γ = 1 :
eta_1_a=thresh_cross_a_eta(POs_cands,γ=1, as=as_LS, η_1_bounds=η_1_bounds_LS,Δt=0.05,timestepUPOs=0.05)
#makePlot_eta_a(eta_1_a,1)
# γ = 5 :
#eta_1_a_gamma5=thresh_cross_a_eta(POs_cands,γ=5, as=as_LS, η_1_bounds=(0.8,1.6),Δt=1.,timestepUPOs=1.,tmax=350)
#makePlot_eta_a(eta_1_a_gamma5,5)
# infinite γ
# The Lorenz-POs run much faster than the Stommel system, and thus, the Stommel system only 
# "sees" the mean of the POs' forcing.
eta_1_a_inf_gamma=thresh_cross_a_eta(POs_cands, po_stats[:,1]; γ=1, as=as_LS, η_1_bounds=η_1_bounds_LS)
makePlot_eta_a(eta_1_a_inf_gamma,"∞")




###################
# non-UPO forcing #
###################
# for each eta and a, use a single new orbit starting at the final point of the previous orbit. 
# Fix the total integration time and then do the plot again with a longer integration time

# Lorenz-Stommel system forced by general Lorenz trajectory
function yellow_blue_plotvals_timeLS(ηs, as, tmax;
    γ=1,
    Δt=0.001,
    t_interm=5
    )

    inits = [0.1,0.1,25.1,1.7,1.0]
    tipped_aη=zeros(length(ηs),length(as))
    for i in 1:length(ηs)
        for j in 1:length(as)

            # run  Lorenz for t_interm Lorenz-timesteps to decorelater subsequent runs
            p=Parameters_fullSyst(γ,3, ηs[i], 0.3, as[j])
            dStommel = ContinuousDynamicalSystem(fullSystem, inits, p, t0=0)
            traj,trajt = trajectory(dStommel, t_interm/γ; Ttr = t_interm-0.01)

            inits = [traj[end][1],traj[end][2],traj[end][3],1.7,1.0]
            dStommel = ContinuousDynamicalSystem(fullSystem, inits, p, t0=0)
            while tipped_aη[i,j]==0.0 && current_time(dStommel) < tmax
                step!(dStommel,Δt,true)
                Ψ = get_state(dStommel)[4]-get_state(dStommel)[5]
                if Ψ < 0.1 
                    tipped_aη[i,j] = current_time(dStommel)
                end
            end
            inits = [get_state(dStommel)[1],get_state(dStommel)[2],get_state(dStommel)[3],1.7,1.0]
        end
    end
    return tipped_aη
end
etavals = 0.8:0.01:1.85
avals_LS = 0.0:0.00125:0.08
LSyb_plotvals_g0_01_tm150=yellow_blue_plotvals_timeLS(0.6:0.01:1.85, avals_LS, 150, γ=0.01)
LSyb_plotvals_g0_05_tm150=yellow_blue_plotvals_timeLS(etavals, avals_LS, 150, γ=0.05)
LSyb_plotvals_g0_1_tm150=yellow_blue_plotvals_timeLS(etavals, avals_LS, 150, γ=0.1)
LSyb_plotvals_g1_tm150=yellow_blue_plotvals_timeLS(etavals, avals_LS, 150, γ=1)

# sort into four time intervals:
function LSsortybplotvals4(yb_vals,upper,mid,lower)
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
LScat_plotvals0_01=LSsortybplotvals4(LSyb_plotvals_g0_01_tm150,150,0.05,0.011)
LScat_plotvals0_05=LSsortybplotvals4(LSyb_plotvals_g0_05_tm150,150,0.22,0.061)
LScat_plotvals0_1=LSsortybplotvals4(LSyb_plotvals_g0_1_tm150,150,0.7,0.125)
LScat_plotvals1=LSsortybplotvals4(LSyb_plotvals_g1_tm150,150,4.0,1.25)


####################################
# combine the heatmap with PO plots
####################################

# extend heatmap values to the size needed
function LSextend_yb_plotvals(yb_plotvals,etavals)
    etastep=round((etavals[end]-etavals[1])/length(etavals),digits=3)
    ext_etavals=0.4:etastep:1.86

    ext_yb_plotvals=zeros(length(avals_LS),length(ext_etavals))
    ext_yb_plotvals[:,1:length(0.4:etastep:etavals[1])].=0.0
    ext_yb_plotvals[:,length(0.4:etastep:etavals[1])+1 : length(0.4:etastep:etavals[1])+length(etavals)].=yb_plotvals'
    ext_yb_plotvals[:,length(0.4:etastep:etavals[1])+1+length(etavals):end].=1.0
    
    return ext_yb_plotvals,ext_etavals
end
ext_hm_vals, ext_etavals = LSextend_yb_plotvals(LScat_plotvals0_01,0.6:0.01:1.85)

function LSmake_overlayed_plot_overlayed0_01(hm_vals, etavals)
    ext_hm_vals = LSextend_yb_plotvals(hm_vals,etavals)

    P1_gamma0_01=plot(xlabel=L"\eta_1", ylabel=L"a",title=L"\gamma = 0.01");
    heatmap!(P1_gamma0_01,ext_etavals, avals_LS, ext_hm_vals,colorbar=false,color=cgrad(:binary,rev=true),alpha=1.);
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
        plot!(P1_gamma0_01,eta_1_a_gamma0_01[pocounter][:,1],eta_1_a_gamma0_01[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            xlim=(0.4,1.85),# (1.028,1.083),  
            ylim=(0.0,0.081),# (0.0325,0.046),
            linewidth=2,
            color=col_c,
            linestyle= line_style_c,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=2,
            xguidefontsize=20,yguidefontsize=20,
            xtickfontsize=12,ytickfontsize=12,
            titlefont=(20,"Computer Modern"),
            legendfontsize=12,leg=:bottomleft)
    end
    #display(P1_gamma0_01)
    return P1_gamma0_01
end
LSmake_overlayed_plot_overlayed0_01(LScat_plotvals0_01,0.6:0.01:1.85)

function LSmake_overlayed_plot_overlayed0_05(hm_vals,etavals)
    ext_hm_vals = LSextend_yb_plotvals(hm_vals,etavals)
    #heatmap(ext_etavals, avals, ext_yb_plotvals_g0_375_tm100,xlabel=L"\eta_1", ylabel=L"a",size=(750,600), title="tmax=50",xlim=(0.4,1.6))

    P1_gamma0_05=plot(xlabel=L"\eta_1", ylabel=L"a",title=L"\gamma = 0.05");
    heatmap!(P1_gamma0_05,ext_etavals, avals_LS, ext_hm_vals,colorbar=false,color=cgrad(:binary,rev=true),alpha=1.);
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
        plot!(P1_gamma0_05,eta_1_a_gamma0_05[pocounter][:,1],eta_1_a_gamma0_05[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            xlim=(0.4,1.85),# (1.028,1.083),  
            ylim=(0.0,0.081),# (0.0325,0.046),
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
    #display(P1_gamma0_375)
    return P1_gamma0_05
end
LSmake_overlayed_plot_overlayed0_05(LScat_plotvals0_05,etavals)

function LSmake_overlayed_plot_overlayed0_1(hm_vals,etavals)
    ext_hm_vals = LSextend_yb_plotvals(hm_vals,etavals)
    #heatmap(ext_etavals, avals, ext_yb_plotvals_g0_375_tm100,xlabel=L"\eta_1", ylabel=L"a",size=(750,600), title="tmax=50",xlim=(0.4,1.6))

    P1_gamma0_1=plot(xlabel=L"\eta_1", ylabel=L"a",title=L"\gamma = 0.1");
    heatmap!(P1_gamma0_1,ext_etavals, avals_LS, ext_hm_vals,colorbar=false,color=cgrad(:binary,rev=true),alpha=1.);
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
        plot!(P1_gamma0_1,eta_1_a_gamma0_1[pocounter][:,1],eta_1_a_gamma0_1[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            xlim=(0.4,1.85),# (1.028,1.083),  
            ylim=(0.0,0.081),# (0.0325,0.046),
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
LSmake_overlayed_plot_overlayed0_1(LScat_plotvals0_1,etavals)

function LSmake_overlayed_plot_overlayed1(hm_vals,etavals)
    ext_hm_vals = LSextend_yb_plotvals(hm_vals,etavals)
    #heatmap(ext_etavals, avals, ext_yb_plotvals_g0_375_tm100,xlabel=L"\eta_1", ylabel=L"a",size=(750,600), title="tmax=50",xlim=(0.4,1.6))

    P1=plot(xlabel=L"\eta_1", ylabel=L"a",title=L"\gamma = 1");
    heatmap!(P1,ext_etavals, avals_LS, ext_hm_vals,colorbar=false,color=cgrad(:binary,rev=true),alpha=1);
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
        plot!(P1,eta_1_a[pocounter][:,1],eta_1_a[pocounter][:,2],
            size=(750,600),
            label=L"PO \ %$pocounter",
            xlim=(0.4,1.85),
            ylim=(0,0.081),
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
    #display(P1_gamma0_375)
    return P1
end
LSmake_overlayed_plot_overlayed1(LScat_plotvals1,etavals)







##########################################################
##########################################################
##########################################################
#                                                        #   
#           dynamics tipping window                      #
#                                                        #
##########################################################
##########################################################
##########################################################


# fix a, start system with small eta and ramp eta to larger values 
# with different speeds, check when the system tips
# repeat this for different random lorenz forcings
mutable struct Parameters5
    γ::Float64
    r::Float64
    ξ_1::Float64
    ζ::Float64
    a::Float64
end

# define a function η(r,t):
η(r,t)=r*t + 0.9

function G4(x, p, t)
    σ = 10.0; ρ = 28.0; β = 8/3;

    dx1 = σ*(x[2] - x[1])
    dx2 = x[1]*(ρ - x[3]) - x[2]
    dx3 = x[1]*x[2] -β*x[3]
    dx4 = ((p.ξ_1+p.a*x[1])-x[4]*(1+abs(x[4]-x[5])))/p.γ
    dx5 = ((η(p.r,t)+p.a*x[1])-x[5]*(p.ζ+abs(x[4]-x[5])))/p.γ

    return SVector{5}(dx1, dx2, dx3, dx4, dx5)
end

function rate_scatter_vals(ens_members, as, r,tmax;
    γ=1,
    Δt=0.01,
    accuracy=5e-3
    )

    inits = [0.1,0.1,25.1,1.7,1.0]
    tipped_aη=zeros(length(as),ens_members)
    for j in 1:length(as)
        #println(j)
        for i in 1:ens_members
            p=Parameters5(γ,r,3, 0.3, as[j])
            dStommel4 = ContinuousDynamicalSystem(G4, inits, p, t0=0)
            while tipped_aη[j,i]==0.0 && current_time(dStommel4) < tmax
                step!(dStommel4,Δt,true)
                Ψ = get_state(dStommel4)[4]-get_state(dStommel4)[5]
                if Ψ < 0.1 
                    #println(j," ",i)
                    tipped_aη[j,i] = η(p.r,current_time(dStommel4)) 
                end
            end
            inits = [get_state(dStommel4)[1],get_state(dStommel4)[2],get_state(dStommel4)[3],1.7,1.0]
        end
    end
    return tipped_aη
end

r_tmax_tubples=[(0.0005,900),(0.06,300),(0.4,100)]

ens_membersRate = 301
avals_rate = 0.0:0.01:0.08
yb_plotvals_g1_rate1=rate_scatter_vals(ens_membersRate, avals_rate, r_tmax_tubples[1][1],r_tmax_tubples[1][2], γ=1)
yb_plotvals_g1_rate2=rate_scatter_vals(ens_membersRate, avals_rate, r_tmax_tubples[2][1],r_tmax_tubples[2][2], γ=1)
yb_plotvals_g1_rate3=rate_scatter_vals(ens_membersRate, avals_rate, r_tmax_tubples[3][1],r_tmax_tubples[3][2], γ=1)

function rate_make_overlayed_plot()

    P1=plot(xlabel=L"\eta_1", ylabel=L"a",title=L"\gamma = 1",legendfontsize=10,size=(750,600),leg=:bottomleft);

    # UPO1 and UPO2 lines
    col_count=1
    for pocounter in 1:2
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
        plot!(P1,eta_1_a[pocounter][:,1],eta_1_a[pocounter][:,2],
            size=(650,500),
            label=L"UPO \ %$pocounter",
            xlim=(0.4,2.2),
            ylim=(0,0.081),
            linewidth=2,
            color=col_c,
            linestyle= line_style_c,
            markershape=markerstyle_c,
            markerstrokecolor=col_c,
            markersize=2,
            xguidefontsize=20,yguidefontsize=20,
            xtickfontsize=12,ytickfontsize=12,
            titlefont=(20,"Computer Modern"))
    end

    # non-UPO rate-dependent forcing
    for i in 1:length(yb_plotvals_g1_rate1[:,1])
        for j in 1:length(yb_plotvals_g1_rate1[1,:])
            if yb_plotvals_g1_rate1[i,j]>0
                plot!(P1,[yb_plotvals_g1_rate1[i,j]],[avals_rate[i]],seriestype=scatter,color=:lime,markerstrokecolor=:lime,alpha=0.025,label=false)
            end
        end
        if i==1
            plot!(P1,[sort(yb_plotvals_g1_rate1[i,:])[floor(Int64,length(yb_plotvals_g1_rate1[1,:])/2)]],[avals_rate[i]],label=L"r=0.0005, \ tmax=900",seriestype=scatter,color=:lime,markerstrokecolor=:white,alpha=1.0)
        else
            plot!(P1,[sort(yb_plotvals_g1_rate1[i,:])[floor(Int64,length(yb_plotvals_g1_rate1[1,:])/2)]],[avals_rate[i]],label=false,seriestype=scatter,color=:lime,markerstrokecolor=:white,alpha=1.0)

        end
    end

    # non-UPO rate-dependent forcing
    for i in 1:length(yb_plotvals_g1_rate2[:,1])
        for j in 1:length(yb_plotvals_g1_rate2[1,:])
            if yb_plotvals_g1_rate2[i,j]>0
                plot!(P1,[yb_plotvals_g1_rate2[i,j]],[avals_rate[i]],seriestype=scatter,color=:orange,markerstrokecolor=:orange,alpha=0.025,label=false)
            end
        end
        if i==1
            plot!(P1,[sort(yb_plotvals_g1_rate2[i,:])[floor(Int64,length(yb_plotvals_g1_rate2[1,:])/2)]],[avals_rate[i]],label=L"r=0.06, \ tmax=300",seriestype=scatter,color=:orange,markerstrokecolor=:white,alpha=1.0)
        else
            plot!(P1,[sort(yb_plotvals_g1_rate2[i,:])[floor(Int64,length(yb_plotvals_g1_rate2[1,:])/2)]],[avals_rate[i]],label=false,seriestype=scatter,color=:orange,markerstrokecolor=:white,alpha=1.0)
        end
    end

    # non-UPO rate-dependent forcing
    for i in 1:length(yb_plotvals_g1_rate3[:,1])
        for j in 1:length(yb_plotvals_g1_rate3[1,:])
            if yb_plotvals_g1_rate2[i,j]>0
                plot!(P1,[yb_plotvals_g1_rate3[i,j]],[avals_rate[i]],seriestype=scatter,color=:blue,markerstrokecolor=:blue,alpha=0.025,label=false)
            end
        end
        if i==1
            plot!(P1,[sort(yb_plotvals_g1_rate3[i,:])[floor(Int64,length(yb_plotvals_g1_rate3[1,:])/2)]],[avals_rate[i]],label=L"r=0.4, \ tmax=100",seriestype=scatter,color=:blue,markerstrokecolor=:white,alpha=1.0)
        else
            plot!(P1,[sort(yb_plotvals_g1_rate3[i,:])[floor(Int64,length(yb_plotvals_g1_rate3[1,:])/2)]],[avals_rate[i]],label=false,seriestype=scatter,color=:blue,markerstrokecolor=:white,alpha=1.0)
        end
    end

    return P1
end
rate_make_overlayed_plot()
