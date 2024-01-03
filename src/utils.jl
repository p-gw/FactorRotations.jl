function zerodiag!(x)
    for i in diagind(x)
        x[i] = 0.0
    end
    return x
end
