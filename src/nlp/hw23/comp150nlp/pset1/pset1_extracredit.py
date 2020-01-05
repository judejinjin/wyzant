# Name: Navraj Narula
# Class: Natural Language Processing
# Assignment 1: Regular Expressions and Finite State Automata
# Date Due: February 12, 2016

# Extra Credit for #6

class FSA:
    """
    define FSA class for NDFA
    """
    def __init__(self, num_states = 0):
        self.num_states = num_states
        self.transitions = {}
        self.final_states = set()
        
    # modfiy add transtions, keep everything else the same.
    
    def add_transition(self, s, currentState, newState):
        if (s, currentState) in self.transitions.keys():
            self.transitions[(s, currentState)] = list(self.transitions[(s, currentState)])
            self.transitions[(s, currentState)].extend(newState)
        else:
            self.transitions[(s, currentState)] = newState
           
    def set_final_state(self, final):
        self.final_states.add(final)
    
    def lookup(self, current, s):  
        if (s, current) in self.transitions:
            return self.transitions[(s, current)]
        else:
            return None
    
    def is_final(self, state):
        return state in self.final_states
        
def createNewState(input_str, currentState, fsa):
    """ helper """
    if currentState is not None:
        index = currentState[0]
        currentNode = currentState[1]
        
        if index != len(input_str):
            destNodes = fsa.lookup(currentNode,input_str[index])
            if destNodes is not None:
                ret = [(index+1, destNode) for destNode in destNodes]
            else:
                ret = None
            return ret
        else:
            return None
    else:
        return None
        
def acceptState(input_str, searchState, fsa):
    """ helper """
    
    if searchState is not None:
        index = searchState[0]
        currentNode = searchState[1]

    if index == len(input_str):
        if fsa.is_final(currentNode):
            return True
        else:
            return False
    else:
        fsa.is_final(currentNode)
        
def NDRecognize(input_str, fsa):
    
    """
    modify DRecognize for NDRecognize;
    implemented by following example in SLP textbook
    data structure used to justify transitions: stack
    """

    index = 0
    initialState = "0"
    agenda = [(index, initialState)]
    currentState = agenda.pop()
    
    while True:
        
        if acceptState(input_str, currentState, FSA):
            return True
        
        else:
            newStateList = createNewState(input_str, currentState, fsa)
        
            if newStateList is not None:
                agenda.extend(newStateList)
                agenda = list(set(agenda))
        
        if not agenda:
            return False
            
        else:
            currentState = agenda.pop()
        
        
def NDRecognizeMulti(input_str, fsa_list): 
    """
    modify DRecognize into NDRecognize, so as to accept dates such as 2/2/2016
    """

    index = 0
    initialState = "0"
    agenda = [(index, initialState)]
    currentState = agenda.pop()
    
    while True:
        
        if acceptState(input_str, currentState, fsa_list[0]):
            if len(fsa_list) == 1:
                return True
            else:
                return NDRecognizeMulti(input_str[currentState[0]:], fsa_list[1:]) or NDRecognizeMulti(input_str[currentState[0]+1:], fsa_list[1:])
        
        else:
            
            newStateList = createNewState(input_str, currentState, fsa_list[0])
            
            if newStateList is not None:
                agenda.extend(newStateList)
                agenda = list(set(agenda))
        
        if not agenda:
            return False
        else:
            currentState = agenda.pop()
                
months = FSA(4)

months.add_transition("0", "0", "1")
months.add_transition("1", "0" ,"2")

num = [str(n) for n in range(1,10)]
map(lambda n: months.add_transition(n, "2", "3"), num)

map(lambda n: months.add_transition(n, "0", "3"), num)

num = [str(n) for n in range(0,3)]
map(lambda n: months.add_transition(n, "2", "3"), num)

months.set_final_state("3")

days = FSA(5)

days.add_transition("0", "0", "1")
days.add_transition("1", "0", "2")
days.add_transition("2", "0", "2")
days.add_transition("3", "0", "3")

num = [str(n) for n in range(1,10)]
map(lambda n: days.add_transition(n, "1", "4"), num)
map(lambda n: days.add_transition(n, "0", "4"), num)

num = [str(n) for n in range(0,10)]
map(lambda n: days.add_transition(n, "2", "4"), num)

days.add_transition("0", "3", "4")
days.add_transition("1", "3", "4")

days.set_final_state("4")

years = FSA("6")

years.add_transition("1", "0", "1")
years.add_transition("2", "0", "2")

years.add_transition("9", "1", "3")
years.add_transition("0", "2", "3")


num = [str(n) for n in [0,1,2,3,4,5,6,7,8,9]]
map(lambda n: years.add_transition(n, "3", "4"),num)

num = [str(n) for n in [0,1,2,3,4,5,6,7,8,9]]
map(lambda n: years.add_transition(n, "4", "5"),num)

years.set_final_state("5")

seps = FSA(2)

seps.add_transition("/", "0", "1")
seps.add_transition(" ", "0", "1")
seps.add_transition("-", "0", "1")

seps.set_final_state("1")

def Test(months, days, years, seps):
    print "\nTest Months FSA"
    for input in ["", "00", "1", "9", "10", "11", "12", "13", "123"]:
        print "'%s'\t%s" %(input, NDRecognizeMulti(input, [months]))
    print "\nTest Days FSA"
    for input in ["", "00", "1", "3", "9", "10", "11", "21", "31", "32", "123"]:
        print "'%s'\t%s" %(input, NDRecognizeMulti(input, [days]))
    print "\nTest Years FSA"
    for input in ["", "1899", "1900", "1901", "1999", "2000", "2001", "2099", "2100"]:
        print "'%s'\t%s" %(input, NDRecognizeMulti(input, [years]))
    print "\nTest Separators FSA"
    for input in ["", ",", " ", "-", "/", "//", ":"]:
        print "'%s'\t%s" %(input, NDRecognizeMulti(input, [seps]))
    print "\nTest Date Expressions FSA"
    for input in ["", "2 31 2000", "3/3/2000", "2-31-2000", "12:31:2000", 
                  "00-31-2000", "12-00-2000", "12-31-0000", 
                  "12-32-1987", "13-31-1987", "12-31-2150"]:
        print "'%s'\t%s" %(input, 
                           NDRecognizeMulti(input, [months, seps, days, seps, years]))

Test(months, days, years, seps)




        
        
        
        
        

    
    

   