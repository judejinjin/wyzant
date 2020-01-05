import ConfigParser
from scapy.all import *

Config = ConfigParser.ConfigParser()
Config.read("inpartb.txt")

src_ip = Config.get("config", "src_ip")
src_port = Config.get("config", "src_port")
dest_ip = Config.get("config", "dest_ip")
dest_port = Config.get("config", "dest_port")
pay_load = Config.get("config", "pay_load")

print "src ip: " + src_ip
print "src port: " + src_port
print "dest ip: " + dest_ip
print "dest port: " + dest_port
print "pay load: " + pay_load

a=IP()
a.src=src_ip
a.dst=dest_ip

b=TCP()
b.sport=int(36000)
b.dport=int(80)

p = a/b/Raw(load=pay_load)

p.show()

print str(p)

print "Sending HTTP packet...\n"

send(p)

print "done."
