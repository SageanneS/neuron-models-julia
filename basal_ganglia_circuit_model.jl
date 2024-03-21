## Backend script for the motor circuit (CBGTC Loop)

using ModelingToolkit, DifferentialEquations, Plots

@parameters t
D = Differential(t)

function NeuralMass(;name, τ=τ, H=H, λ=λ, r=r)
       
       sts    = @variables x(t)=1.0 y(t)=1.0 jcn(t)=0.0
       params = @parameters τ=τ H=H λ=λ r=r

       eqs = [D(x) ~ y - ((2/τ)*x),
              D(y) ~ -x/(τ*τ) + (H/τ)*((2*λ)/(1 + exp(-r*(jcn))) - λ)]
       
       return ODESystem(eqs, t, sts, params; name=name)
end

function Connections(;name, sys=sys)
       
       params =  @parameters C_Cor=60 C_BG_Th=60 C_Cor_BG_Th=5 C_BG_Th_Cor=5

       adj_matrix = [0 0 0 0 0 0 0 0; 
       -0.5*C_BG_Th*sys[1].x -0.5*C_BG_Th*sys[2].x C_BG_Th*sys[3].x 0 0 0 0 0;
                     0 -0.5*C_BG_Th*sys[2].x 0 0 0 0 C_Cor_BG_Th*sys[7].x 0;
                     0 -0.5*C_BG_Th*sys[2].x C_BG_Th*sys[3].x 0 0 0 0 0;
                     0 0 0 -0.5*C_BG_Th*sys[4].x 0 0 0 0;
                     0 0 0 0 C_BG_Th_Cor*sys[5].x 0 6*C_Cor*sys[7].x 0;
                     0 0 0 0 0 4.8*C_Cor*sys[6].x 0 -1.5*C_Cor*sys[8].x;
                     0 0 0 0 0 0 1.5*C_Cor*sys[7].x 3.3*C_Cor*sys[8].x]
        begin
               eqs = []
               for region_num in 1:length(sys)
                      push!(eqs, sys[region_num].jcn ~ sum(adj_matrix[region_num,:]))
               end
        end

        return @named Circuit = ODESystem(eqs, systems = sys)
end

function CreateCircuit(;name)
       
       # Create Regions
       @named Str = NeuralMass(τ=0.0022, H=20, λ=300, r=0.3)
       @named GPe = NeuralMass(τ=0.04, H=20, λ=400, r=0.1)
       @named STN = NeuralMass(τ=0.01, H=20, λ=500, r=0.1)
       @named GPi = NeuralMass(τ=0.014, H=20, λ=400, r=0.1)
       @named Th  = NeuralMass(τ=0.002, H=10, λ=20, r=5)
       @named EI  = NeuralMass(τ=0.01, H=20, λ=5, r=5)
       @named PY  = NeuralMass(τ=0.001, H=20, λ=5, r=0.15)
       @named II  = NeuralMass(τ=2.0, H=60, λ=5, r=5)

       # Connect Regions
       sys = [Str, GPe, STN, GPi, Th, EI, PY, II]
       @named Circuit = Connections(sys=sys)

end
@named CBGTC_Circuit = CreateCircuit()

sim_dur = 10.0 # Simulate for 10 Seconds
prob = ODAEProblem(structural_simplify(CBGTC_Circuit), [], (0.0, sim_dur), [])
sol = solve(prob)

## Plot the Solution:
p1 = plot(sol.t, sol[3,:], lc=:purple, xlabel = "time (s)", ylabel="GPe.x", label="GPe.x(t)")
p2 = plot(sol[3,:], sol[4,:], lc =:blue, xlabel = "GPe.x", ylabel="GPe.y", label = "GPe.x(t), GPe.y(t)", title = "GPe Region", titlefont=12)
p3 = plot(sol.t, sol[5,:], lc=:cyan, xlabel = "time (s)", ylabel="STN.x", label="STN.x(t)")
p4 = plot(sol[5,:], sol[6,:], lc =:blue, xlabel = "STN.x", ylabel="STN.y", label = "STN.x(t), STN.y(t)", title = "STN Region", titlefont=12)
p5 = plot(sol.t, sol[7,:], lc=:red, xlabel = "time (s)", ylabel="GPi.x", label="GPi.x(t)")
p6 = plot(sol[7,:], sol[8,:], lc =:blue, xlabel = "GPi.x", ylabel="GPi.y", label = "GPi.x(t), GPi.y(t)", title = "GPi Region", titlefont=12)
p7 = plot(sol.t, sol[9,:], lc=:green, xlabel = "time (s)", ylabel="Th.x", label="Th.x(t)", ylim=(-1,1))
p8 = plot(sol[9,:], sol[10,:], lc =:blue, xlabel = "Th.x", ylabel="Th.y", label = "Th.x(t), Th.y(t)", title = "Thalamus Region", titlefont=12)
p9 = plot(sol.t, sol[13,:], lc=:orange, xlabel = "time (s)", ylabel="PY.x", label="PY.x(t)")
p10 = plot(sol[13,:], sol[14,:], lc =:blue, xlabel = "PY.x", ylabel="PY.y", label = "PY.x(t), PY.y(t)", title = "Pyramidal Region", titlefont=12)

l_oscillations = grid(5, 2, heights=[0.2, 0.2 ,0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2])
plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, layout = l_oscillations, size=(800,800), lw=2.5)
