# @kris96tian

 
class FmIndex:
    def __init__(self, text):
        self.text = text
        self.suffix_array = self.build_suffix_array()
        self.bwt = self.construct_bwt()
        self.first_occurrence, self.count = self.create_fm_index()

    def build_suffix_array(self):
        suffixes = []
        for i in range(len(self.text)):
            suffixes.append((self.text[i:], i))
        suffixes.sort()
        suffix_array = []
        for suffix in suffixes:
            suffix_array.append(suffix[1])
        return suffix_array

    def construct_bwt(self):
        bwt = ''
        for i in self.suffix_array:
            if i > 0:
                bwt += self.text[i - 1]
            else:
                bwt += self.text[-1]
        return bwt


    def create_fm_index(self):
        first_occurrence = {}
        count = {}
        symbols = set(self.bwt)
        sorted_bwt = sorted(self.bwt)
        for symbol in symbols:
            first_occurrence[symbol] = sorted_bwt.index(symbol)   
        for symbol in symbols:
            count[symbol] = [0] * (len(self.bwt) + 1)    
        for i, symbol in enumerate(self.bwt): # Calc. count for each symbol
            for s in symbols:
                count[s][i + 1] = count[s][i] + 1 if s == symbol else count[s][i]
        return first_occurrence, count


    def fm_substring_search(self, pattern):
        top = 0
        bottom = len(self.bwt) - 1
        for i in range(len(pattern) - 1, -1, -1):
            symbol = pattern[i]
            top = self.first_occurrence[symbol] + self.count[symbol][top]
            bottom = self.first_occurrence[symbol] + self.count[symbol][bottom + 1] - 1
            if top > bottom:
                return None
        return (top, bottom)

    def display_fm_index(self):
        print("Suffix Array:", self.suffix_array)
        print("BWT:", self.bwt)
        print("First Occurrence:", self.first_occurrence)
        print("Count:")
        for symbol, counts in self.count.items():
            print(f"    {symbol}: {counts}")
