'''
string modification
'''

def get_char_series_by_two(chars):
    '''
    chars should have two character
    '''
    ch1 = chars[0]
    ch2 = chars[-1]
    lchar = [chr(ch) for ch in range(ord(ch1), ord(ch2)+1)]
    return lchar

