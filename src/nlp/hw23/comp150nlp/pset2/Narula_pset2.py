import sys
from collections import defaultdict
from math import log, exp
from nltk.corpus import brown
import numpy as np

unknown_token = "<UNK>"  # unknown word token.
start_token = "<S>"  # sentence boundary token.
end_token = "</S>"  # sentence boundary token.

""" Implement any helper functions here, e.g., for text preprocessing.
"""

class BigramLM:
    def __init__(self, vocabulary = set()):
        
        self.vocabulary = vocabulary
        self.unigram_counts = defaultdict(float)
        self.bigram_counts = defaultdict(float)
        self.log_probs = {}
        
    """ Implement the functions EstimateBigrams, CheckDistribution, Perplexity and any 
    other function you might need here.
    """
        
    def count_unigrams(self, sentences):      
        self.unigram_counts = count_occurences(sentences, self.vocabulary)
            
    def count_bigrams(self, sentences):
            
        bigramm_counts = defaultdict(float)
        for sentence in sentences:
            for index, word in enumerate(sentence):
                if index > 0:
                    w1, w2 = sentence[index-1], sentence[index]
                    if w1 not in self.vocabulary:
                        w1 = unknown_token
                    if w2 not in self.vocabulary:
                        w2 = unknown_token #tokens not in vocab, unknown token
                    self.bigram_counts[(w1,w2)] += 1
        return bigramm_counts
                                                
    def logprob(self):
        for bigram, bigramcount in self.bigram_counts.iteritems():
            u, v = bigram
            if self.unigram_counts[u] == 0:
                print u
            self.log_probs[bigram] = log(bigramcount/self.unigram_counts[u])
                

    def EstimateBigrams(self, sentences):
        BigramLM.count_unigrams(self,sentences)
        BigramLM.count_bigrams(self,sentences)
        BigramLM.logprob(self)
            
                   
    def CheckDistribution(self):     
        prob_sum = defaultdict(float)
            
        for key, val in self.log_probs.iteritems():
            prob_sum[key[0]] += exp(val)
                
        for key, val in prob_sum.iteritems():
            assert np.isclose(val, 1.0)
                
    def Perplexity(self, sentences):
        
        logp = 0.0
        numwords = 0.0
            
        for sentence in sentences:
            for index, word in enumerate(sentence):
                numwords += 1
                if index > 0:
                    w1, w2 = sentence[index-1], sentence[index]
                    if w1 not in self.vocabulary:
                        w1 = unknown_token
                    if w2 not in self.vocabulary:
                        w2 = unknown_token
                    
                    key = (w1,w2)
                    
#                    if key not in self.log_probs.keys():
#                        # takes care of key errors <<replace>>
#                        logp += 0#
#                        #self.log_probs[(unknown_token, unknown_token)]
#                        
#                    else:
#                        logp += self.log_probs[key]
                    
                    try:
                        logp += self.log_probs[key]
                    except KeyError:
                        logp += 0 
                
        logp *= -1/numwords
        return exp(logp)
        
            
    def laplace_smoothing(self, sentences):
        voc_size = len(self.vocabulary)
        for bigram, bigramcount in self.bigram_counts.iteritems():
            u, v = bigram
            self.log_probs[bigram] = log((bigramcount + 1.0)/(self.unigram_counts[u] + voc_size))
            
    def simple_linear_interpolation(self, sentences):
        
        l1, l2 = [0.5] * 2
        
        wordcount = count_words(sentences)
        
        for key, val in self.bigram_counts.iteritems():
            w,v = key
            bigram_prob = val / self.unigram_counts[w]
            unigram_prob = self.unigram_counts[v] / wordcount
            self.log_probs[key] = log(l2 * bigram_prob + l1 * unigram_prob)
            
    def deleted_interpolation(self, held_out_sentences, sentences):
        
        l1, l2 = [0] * 2
        held_out_count = count_words(held_out_sentences)
        
        
        bi_held_counts = self.count_bigrams(held_out_sentences)
        uni_held_counts = self.count_unigrams(held_out_sentences)
        
        for bigram, bicount in self.bigram_counts.iteritems():
            w,v = bigram
            
            
            if bi_held_counts[bigram] > 0:
                 print bi_held_counts
                
                
                 bitest = (bi_held_counts[bigram] - 1) / (uni_held_counts[w] - 1) if (uni_held_counts[w] - 1) != 0 else 0
                 unitest = (uni_held_counts[v] - 1) / (held_out_count - 1) if (held_out_count - 1) != 0 else 0
                
                 if bitest > unitest:
                     l2 += bi_held_counts[bigram]
                 else:
                     l1 += bi_held_counts[bigram]
                     
     
        l1_norm = l1 / (l1 + l2) if l1 + l2 != 0 else 0
        l2_norm = l2 / (l1 + l2) if l1 + l2 != 0 else 0
        
        print "The estimated weights are: ", l1_norm, " and ", l2_norm

def count_words(sentences):
        
    word_count = 0 
    for sentence in sentences:
        word_count += len(sentence)

    return word_count       
                    
    
def create_vocabulary(sentences):
    
    voc = set()
    for sentence in sentences:
        for word in sentence: #token
            voc.add(word)
    voc.add(unknown_token)
    voc.add(start_token)
    voc.add(end_token)
    return voc
        
def count_occurences(sentences, vocabulary):
    
    word_counts = defaultdict(float)
    for sentence in sentences:
        for word in sentence:
            if word not in vocabulary:
                word_counts[unknown_token] += 1
            else:
                word_counts[word] += 1
    return word_counts
    
      
def PreprocessText(sentences, vocabulary):
    
    new_sentences = sentences
    
    word_counts = count_occurences(new_sentences, vocabulary)
    
    for new_sentence in new_sentences:
        for index, word in enumerate(new_sentence): 
            
            if word_counts[word] < 2 or word not in vocabulary:
                new_sentence[index] = unknown_token
                            
    new_sentences = [[start_token] + sentence + [end_token] for sentence in new_sentences]
    
            
    return new_sentences
    
   
def main():
    training_set = brown.sents()[:50000]
    held_out_set = brown.sents()[-6000:-3000]
    test_set = brown.sents()[-3000:]

    vocabulary = create_vocabulary(training_set)
   
    """ Transform the data sets by eliminating unknown words and adding sentence boundary 
    tokens.
    """
    
    training_set_prep = PreprocessText(training_set, vocabulary)
    held_out_set_prep = PreprocessText(held_out_set, vocabulary)
    test_set_prep = PreprocessText(test_set, vocabulary)

    """ Print the first sentence of each data set.
    """
    print training_set_prep[0]
    print
    print held_out_set_prep[0]
    print
    print test_set_prep[0]
    print

    """ Estimate a bigram_lm object, check its distribution, compute its perplexity.
    """
    
    blm = BigramLM(vocabulary)
    blm.EstimateBigrams(training_set_prep)
    blm.CheckDistribution()
    print "The perplexity of the training set is: ", blm.Perplexity(test_set_prep)
    
    

    """ Print out perplexity after Laplace smoothing.
    
    """ 
    blm.EstimateBigrams(training_set_prep)
    blm.laplace_smoothing(training_set_prep)
    print "The perplexity of the training set after Laplace smoothing is: ", blm.Perplexity(test_set_prep)
    

    """ Print out perplexity after simple linear interpolation (SLI) with lambda = 0.5.
    """ 
    
    blm.EstimateBigrams(training_set_prep)
    blm.simple_linear_interpolation(training_set_prep)
    print "The perplexity of the training set after SLI is: ", blm.Perplexity(test_set_prep)
    

    """ Estimate interpolation weights using the deleted interpolation algorithm on the 
    held out set and print out. #held-out 
    """ 
    
    blm.EstimateBigrams(training_set_prep)
    blm.deleted_interpolation(held_out_set_prep, training_set_prep)

    """ Print out perplexity after simple linear interpolation (SLI) with the estimated
    interpolation weights.
    """ 
    
    print "The perplexity of the training set with estimated weights is: ", blm.Perplexity(test_set_prep)

if __name__ == "__main__": 
    main()







    