using TimerOutputs
using Printf

global GLOBAL_TIMER = TimerOutput()
global ENABLE_PROFILING = false
    
macro trace(name, expr)
    esc(quote
        if ENABLE_PROFILING
            @timeit GLOBAL_TIMER $name begin
                $expr  # 🔹 Ensure the expression runs in the correct scope
            end
        else
            $expr  # 🔹 Run normally without profiling
        end
    end)
end


function print_timings(TO :: TimerOutput)
    if ENABLE_PROFILING
        println("=== Performance Metrics ===")
        show(TO; allocations=false)  # Print timing results
    end
end
