import sys
from scapy.all import *

if len(sys.argv) != 2:
    print "usage: qa.py ip_address\n e.g. qa.py 10.0.0.9"
    sys.exit(1)

p = IP(dst=sys.argv[1])/ICMP()

print "Sending ICMP packet...\n"

(ans,unans) = sr(p)

print "\nReceived " + str(len(ans)) + " replies:\n"

if len(ans) > 0:
    print ans

if len(unans) > 0:
    print unans
