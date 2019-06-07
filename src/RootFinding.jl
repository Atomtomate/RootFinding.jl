module RootFinding
# Source from mmas.github.io
export brent, rf

function brent(f::Function, xInit::T1,  Δx = 2.0, xtol::T2 = 2eps(T1), ytol::T2 = 2eps(T1), maxiter=100) where {T1 <: Real, T2 <: Real}
    x0 = xInit - Δx
    x1 = xInit + Δx
    while maxiter > 0
        y0 = f(x0)
        y1 = f(x1)
        if sign(y0) != sign(y1)     # Root is bracketed, increade search spcae
            break
        end
        maxiter -= 1
        x0 -= Δx
        x1 += Δx
    end
    if abs(y0) < abs(y1)
        x0, x1 = x1, x0
        y0, y1 = y1, y0
    end
    x2 = x0
    y2 = y0
    x3 = x2
    use_bisection = true
    machine_eps = eps(Float64)

    for i in 1:maxiter
        if abs(x1 - x0) < xtol
            print(i)
            print(" (brent): ")
            println(x1)
            return x1
        end

        # inverse quad. interp., otherwise secant
        if abs(y0 - y2) > ytol && abs(y1 - y2) > ytol
            x = x0*y1*y2 / ((y0-y1)*(y0-y2)) + 
                x1*y0*y2 / ((y1-y0)*(y1-y2)) +
                x2*y0*y1 / ((y2-y0)*(y2-y1))
        else
            x = x1 - y1*(x1-x0)/(y1-y0)
        end

        # bisection
        δ = abs(machine_eps*2*abs(x1))
        min1 = abs(x-x1)
        min2 = abs(x1-x2)
        min3 = abs(x2-x3)
        if (x < (3.0*x0 + x1)/4.0 && x > x1) || 
                (use_bisection && min1 >= 0.5*min2) || (use_bisection && min2 < δ) ||
                (!use_bisection && min1 >= 0.5*min3) || (!use_bisection && min3 < δ)
            x = 0.5*(x0+x1)
            use_bisection = true
        else
            use_bisection = false
        end

        y = f(x)
        if abs(y) < ytol
            print(i)
            print(" (brent): ")
            println(x)
            return x
        end
        x3 = x2
        x2 = x1

        if sign(y0) != sign(y)
            x1 = x
            y1 = y
        else
            x0 = x
            y0 = y
        end
        if abs(y0) < abs(y1)
            x0, x1 = x1, x0
            y0, y1 = y1, y0
        end
    end
    throw(InterruptException("\nBrent's Method did not converge after given number of iterations"))
end

function rf(f::Function, x0::T1, x1::T1, xtol::T2 = 2eps(T1), ytol::T2 = 2eps(T1), maxiter=1000) where {T1 <: Real, T2 <: Real}
    y0 = f(x0)
    y1 = f(x1)

    for i in 1:maxiter
        x = x1 - y1*(x1-x0)/(y1-y0)
        if min(abs(x-x0), abs(x-x1)) < xtol
            print(i)
            print(" (rf): ")
            println(x)
            return x
        end
        y = f(x)
        if abs(y) < ytol
            print(i)
            print(" (rf): ")
            println(x)
            return x
        end
        if sign(y0*y) == 1
            x0 = x
            y0 = y
        else
            x1 = x
            y1 = y
        end
    end
end

end # module
