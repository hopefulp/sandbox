### image extraction from total images

from my_arith import divide_int
import sys


class Images:
    """
    get total images,
    divide total images into train-test sets ('t-t' type)
    or train-validation-test sets ('t-v-t' type)
    """
    def __init__(self, total_images, nsets=None):
        self.total_images = total_images
        self.nsets=nsets
        ### in case nsets is defined
        if isinstance(nsets, int):
            if self.nsets == 1:
                self.train_all=self.total_images[:]
            else:
                ### defin self.parts
                self.parts = divide_int(len(self.total_images), self.nsets)

    # these method works independently
    def get_training_images(self, d_list=None):
        print(d_list)
        if d_list:
            return self.total_images[d_list[0]:d_list[1]]
        elif self.nsets == 1:
            return self.train_all
        else:
            print(self.parts, len(self.total_images))
            return self.total_images[:self.parts[-2]]

    def get_test_images(self, d_list=None):
        if d_list:
            return self.total_images[d_list[0]:d_list[1]]
        elif self.nsets == 1:
            print("error: the whole set cann't be test sets")
            sys.exit(10)
        else:
            print(self.parts, len(self.total_images))
            return self.total_images[self.parts[-2]:]


class Images1:
    """
    get total images,
    divide total images into train-test sets ('t-t' type)
    or train-validation-test sets ('t-v-t' type)
    """
    def __init__(self, total_images, nsets, job):
        self.total_images = total_images
        self.nsets=nsets
        if self.nsets == 1:
            self.train_all=self.total_images[:]
        else:
            self.parts = divide_int(len(self.total_images), self.nsets)

    # these method works independently
    def get_training_images(self):
        if self.nsets == 1:
            return self.train_all
        else:
            print(self.parts, len(self.total_images))
            return self.total_images[:self.parts[-2]]

    def get_test_images(self):
        if self.nsets == 1:
            print("error: the whole set cann't be test sets")
            sys.exit(10)
        else:
            print(self.parts, len(self.total_images))
            return self.total_images[self.parts[-2]:]

    def get_val_train_images(self, n):
        parts = divide_int(len(self.total_images), self.nsets)
        print(parts, len(self.total_images))
        if n == 0:
            v_start = 0         # v for valence set area
            v_stop  = parts[0]
        else:
            v_start = parts[n-1]
            v_stop  = parts[n]
        val_images = self.total_images[v_start:v_stop]
        if n == 0:
            tra_images = self.total_images[v_stop:parts[-2]]
        else:
            tra_images = self.total_images[:v_start]
            if n != self.nsets-2:
                tra_images.extend(self.total_images[v_stop:parts[-2]])
        return tra_images, val_images

class Images_old:
    """
    get total images,
    divide total images into train-test sets ('t-t' type)
    or train-validation-test sets ('t-v-t' type)
    """
    def __init__(self, total_images):
        self.total_images = total_images
    # these method works independently
    def get_training_images(self, ndiv):
        self.ndiv = ndiv
        self.parts = divide_int(len(self.total_images), self.ndiv)
        print(self.parts, len(self.total_images))
        return self.total_images[:self.parts[-2]]

    def get_test_images(self, ndiv):
        self.ndiv = ndiv
        self.parts = divide_int(len(self.total_images), self.ndiv)
        print(self.parts, len(self.total_images))
        return self.total_images[self.parts[-2]:]

    def get_val_train_images(self, ndiv, n):
        #self.ndiv = ndiv
        parts = divide_int(len(self.total_images), ndiv)
        print(parts, len(self.total_images))
        if n == 0:
            v_start = 0
            v_stop  = parts[0]
        else:
            v_start = parts[n-1]
            v_stop  = parts[n]
        val_images = self.total_images[v_start:v_stop]
        if n == 0:
            tra_images = self.total_images[v_stop:parts[-2]]
        else:
            tra_images = self.total_images[:v_start]
            if n != ndiv-2:
                tra_images.extend(self.total_images[v_stop:parts[-2]])
        return tra_images, val_images
        
    #def         


'''        
class Images_old:
    
    get total images,
    divide total images into train-test sets ('t-t' type)
    or train-validation-test sets ('t-v-t' type)
    
    def __init__(self, total_images, total_ndiv=5):
        self.total_images = total_images
        self.total_ndiv=5
        self.parts = divide_int(len(self.total_images), self.total_ndiv)
        print(self.parts, len(self.total_images))
'''        
