import numpy as np
from common import whereami

def outcar_trim(total_images):
    print(f"Before data process: num_total {len(total_images)} in {whereami()}")
    
    ### remove high values
    new_images = [ image for image in total_images if image.get_potential_energy() < -950 ]
    
    ### remove high stepping energy

    to_delete=[]
    for i, image in enumerate(new_images):
        pot = image.get_potential_energy()
        if i == 0:
            ptmp = pot
        else:
            if pot > ptmp + 0.8:
                to_delete.append(i)
            else:
                ptmp = pot
    to_delete.reverse()     # in line reverse
    print(to_delete)
    for i in to_delete:
        del new_images[i]
        
    print(f"After data process: num_total {len(new_images)} in {whereami()}")
    return new_images
