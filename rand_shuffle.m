% Randomly shuffles the entries in a vector, keeping the size the same.
function Vout = rand_shuffle(Vin)
Vout = Vin(randperm(length(Vin)));
end
