# files = ['/Users/swang/Simulation/potentials/EAM/Ram_song_potA.fs']
files = ['/Users/swang/Simulation/make-potential/ram_potb/Ram_potB.txt']
# files = [ "/Users/swang/Simulation/make-potential/Mendelev.fs" ]
# files = ['/Users/swang/Simulation/potentials/EAM/Ram_potB.fs']

pair = ["* * " + files[0] + " Fe H "]
parameters = {"pair_style": "eam/fs", "pair_coeff": pair}

# def pick_elements(potpara, elems):
#     for e in elems:
#         potpara.parameters["pair_coeff"][0]  += " " + e
