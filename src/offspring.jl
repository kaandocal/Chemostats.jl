"""
    ExactPartition(f = nothing)

Divides `x` among all offspring, such that offspring `i` receives `f[i] * x`. The fractions `f`
must be nonnegative and sum to 1. By default, `x` is divided equally. 
"""
struct ExactPartition{F}
    f::F

    function ExactPartition(f = nothing)
        check_fractions(f)
        new{typeof(f)}(f)
    end
end

"""
    BinomialPartition(f = nothing)

Divides a count `x` (rounded to an integer) among the offspring of a cell, following 
a multinomial distribution. Each molecule is independently assigned to offspring `i` 
with probability `f[i]`. The fractions `f` must be nonnegative and sum to 1. 
By default, all offspring are equally likely.
"""
struct BinomialPartition{F}
    f::F

    function BinomialPartition(f = nothing)
        check_fractions(f)
        new{typeof(f)}(f)
    end
end

"""
    AllOrNothingPartition(f = nothing)

Assigns all of a quantity `x` to a single offspring, chosen at random with probabilities `f`.
The probabilities `f` must be nonnegative and sum to 1. By default, all offspring are
equally likely.
"""
struct AllOrNothingPartition{F}
    f::F

    function AllOrNothingPartition(f = nothing)
        check_fractions(f)
        new{typeof(f)}(f)
    end
end

const Partition = Union{ExactPartition,BinomialPartition,AllOrNothingPartition}

check_fractions(::Nothing) = nothing
function check_fractions(f)
    @argcheck all(>=(0), f) "Fractions must be nonnegative"
    @argcheck sum(f) ≈ 1 "Fractions must sum to 1"
end

"""
    Offspring([n,] rules...)

Describes initial conditions and parameters for the `n` offspring of a dividing cell.
Each rule is of the form `var => value`, where `value` can be one of the following:

* [`ExactPartition(f)`](@ref ExactPartition): `var` is divided exactly among the offspring
* [`BinomialPartition(f)`](@ref BinomialPartition): `var` is partitioned binomially among the offspring
* [`AllOrNothingPartition(f)`](@ref AllOrNothingPartition): `var` goes entirely to one offspring
* `[ val1, val2, ..., valn ]`: the `i`-th offspring gets the `i`-th value
* `val`: each offspring gets `val`,
* `dist::Distribution`: each offspring gets an independent random sample from `dist`

All other variables are copied from the parent. 

An `Offspring` can be used as the `divide` keyword of a [`DECell`](@ref).
Alternatively, `divide` can be a function `int -> offspring` that returns an `Offspring`
object.
"""
struct Offspring{R<:Tuple}
    n::Int
    rules::R

    function Offspring(n::Integer, rules::Pair...)
        @argcheck n >= 0 "Number of offspring must be nonnegative"
        for (var, rule) in rules
            len = rule_length(rule)
            @argcheck isnothing(len) || len == n "Rule for $var has $len entries, but there are $n offspring"
        end
        @argcheck all(var -> var isa Integer, first.(rules)) || !any(var -> var isa Integer, first.(rules)) "Variables must be all symbols or all indices"
        new{typeof(rules)}(n, rules)
    end
end

function Offspring(rules::Pair...)
    lengths = filter(!isnothing, [ rule_length(last(rule)) for rule in rules ])
    Offspring(isempty(lengths) ? 2 : first(lengths), rules...)
end

# Returns itself, so that an `Offspring` can be used as the `divide` function of a cell
(offspring::Offspring)(int) = offspring

# Number of offspring a rule specifies, or `nothing` if it works for any number
rule_length(rule::Partition) = isnothing(rule.f) ? nothing : length(rule.f)
rule_length(rule::AbstractVector) = length(rule)
rule_length(rule) = nothing

fractions(f, n) = isnothing(f) ? fill(1 / n, n) : f

# Values of a variable with value `x` in the parent for each of the `n` offspring
rule_values(rule::ExactPartition, x, n, rng) = fractions(rule.f, n) .* x
rule_values(rule::BinomialPartition, x, n, rng) = multinomial(rng, round(Int, x), fractions(rule.f, n))

function rule_values(rule::AllOrNothingPartition, x, n, rng)
    winner = findfirst(==(1), multinomial(rng, 1, fractions(rule.f, n)))
    [ i == winner ? x : zero(x) for i in 1:n ]
end
rule_values(rule::AbstractVector, x, n, rng) = rule
rule_values(rule::Number, x, n, rng) = fill(rule, n)
rule_values(rule, x, n, rng) = [ rand(rng, rule) for _ in 1:n ]

# Assigns each of `N` items to one of the categories with probabilities `f`
function multinomial(rng, N::Integer, f)
    @argcheck N >= 0 "Cannot partition a negative count"
    counts = zeros(Int, length(f))
    cdf = cumsum(f)
    for _ in 1:N
        i = min(searchsortedfirst(cdf, rand(rng) * cdf[end]), length(f))
        counts[i] += 1
    end
    counts
end

parent_value(int, var::Integer) = int.u[var]
parent_value(int, var) = int[var]

# Converts an `Offspring` into a vector of `(u0, p)` daughters, as returned by the `divide` function of a `DECell`
function offspring_states(offspring::Offspring, int, rng)
    (; n, rules) = offspring
    values = map(rules) do (var, rule)
        x = rule isa Partition ? parent_value(int, var) : nothing
        var => rule_values(rule, x, n, rng)
    end

    map(1:n) do i
        u0 = if all(rule -> first(rule) isa Integer, rules)
            u = copy(int.u)
            for (var, vals) in values
                u[var] = vals[i]
            end
            u
        else
            [ var => vals[i] for (var, vals) in values ]
        end
        (u0 = u0, p = nothing)
    end
end
