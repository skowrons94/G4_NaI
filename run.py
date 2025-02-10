import os
from tqdm import tqdm

# Check if plot direcoty exists, if not create it
if not os.path.exists('output'):
    os.makedirs('output')

# Clear output folder
os.system('rm output/* >/dev/null 2>&1')

sources = ["40K.in", "210Pb.in", "222Ra.in", "228Th.in", "232Th.in", "238U.in"]

for source in tqdm(sources):
    
    # Run the simulation (ignore the output)
    os.system('./rdecay01 mac/{} >/dev/null 2>&1'.format(source))

    # Move the file that contains "_t0.root" to the output folder
    os.system('mv {}.root output/{}.root'.format(source.split('.')[0], source.split('.')[0]))