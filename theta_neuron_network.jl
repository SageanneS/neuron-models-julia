using ModelingToolkit, DifferentialEquations, Plots

@parameters t
D=Differential(t)
function ThetaNeuron(;name, η=η, α_inv=α_inv, k=k)

    params = @parameters η=η α_inv=α_inv k=k
    # η:     Constant drive
    # α_inv: Time to peak of spike
    # k:     All-to-all coupling strength

    sts    = @variables θ(t)=0.0 g(t)=0.0 jcn(t)=0.0
    # θ(t):   Theta neuron state
    # g(t):   Synaptic current
    # jcn(t): Synaptic input

   eqs = [D(θ) ~ 1-cos(θ) + (1+cos(θ))*(η + k*g),
          D(g) ~ α_inv*(jcn - g)]

    return ODESystem(eqs, t, sts, params; name=name)

end

function ThetaNetwork(;name, N=N, η0=η0, Δ=Δ, n=n)

    # N:     Neuron population size
    # η0, Δ: Distribution parameters for generating a constant drive into each neuron
    # n:     Spike shape parameter

    # Create Neuron Network:
    network = [] 
    for i = 1:N
        η  = rand(Cauchy(η0, Δ)) # Constant Drive
        @named neuron = ThetaNeuron(name=Symbol("neuron$i"), η=η, α_inv=1.0, k=-2.0)
        push!(network, neuron)
    end

    a_n = 2.0^n*(factorial(n)^2.0)/(factorial(2*n))
    sysx = [a_n*(1-cos(neuron.θ))^n for neuron in network]
    adj_matrix = ones(N,N)
    adjx = adj_matrix .* sysx

    eqs = []
    for i = 1:N
        push!(eqs, network[i].jcn ~ sum(adjx[i])/N)
    end
    return @named Network = ODESystem(eqs, systems = network)

end
@named ThetaCircuit = ThetaNetwork(N=500, η0=1.0, Δ=0.05, n=3)

sim_dur = 50.0 # Simulate for 10 Seconds
prob = ODAEProblem(structural_simplify(ThetaCircuit), [], (0.0, sim_dur), [])
sol = solve(prob, Tsit5())