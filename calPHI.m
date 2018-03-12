function PHI = calPHI(Q)
    PHI = 0;
    for q=1:Q
        PHI = PHI + calEuler(q);
    end
end