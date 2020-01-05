import sys
import string
import random

from scapy.all import *

def randomTenLetters():
    s = ""
    for n in range(0, 10):
        s += random.choice(string.letters)
    return s
    
if len(sys.argv) != 3:
    print "usage: qc.py src_port dest_port\n e.g. qa.py 36000 3000 "
    sys.exit(1)

src_port = int(sys.argv[1])
dest_port = int(sys.argv[2])


a = IP(dst='127.0.0.1')

b = TCP()
b.sport = src_port



for n in range(0, 21):
    dport = dest_port + n
    print 'sending to dest_port: ' + str(dport)
    b.dport = dport
    p = a/b/""
    send(p)


for n in range(0, 5):
    b.dport = dport
    pay_load = randomTenLetters()
    print 'sending random 10 letters as pay load:' + str(pay_load)
    p = a/b/Raw(load=pay_load)
    send(p)
    
