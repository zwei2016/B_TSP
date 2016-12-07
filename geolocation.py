
import sys
from geopy.geocoders import Nominatim


if __name__=="__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    fin = open(infile, 'r')
    fout = open(outfile, 'w')
    geolocator = Nominatim()
    
    for sent in fin: 
        #name = sent.split('\n')
        location = geolocator.geocode(sent.strip(),timeout=10)
        print(sent.strip(), location.latitude, location.longitude)
        fout.write("%s\t" % sent.strip())
        fout.write( "%s\t" % location.latitude) 
        fout.write( "%s\n" % location.longitude)
    
    fout.close()
    fin.close()