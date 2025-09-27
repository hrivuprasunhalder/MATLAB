function BE= binding_energy(Z,A,M,isAtomicMass)
    mp=1.00727646688;
    mn=1.00866491588;
    me=0.000548579909;
    if isAtomicMass
       M=M-Z*me;
    end
    deltaM=Z*mp+(A-Z)*mn-M;
    BE=deltaM*931.5;
end    