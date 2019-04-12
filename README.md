# r.stream.power - A GRASS GIS module for stream power calculatiion
*r.stream.power* is a GRASS GIS python script able co calculate teh stream power of a river cross section. The script starts from calculating the food discharge, according to the flood discharge handbook provided by River Tiber Authority.

## Run example
./r.discharge_abt.py dem=elevation  time=200 outlets=334736.07613,4747333.21315

### Run example with r.stream.power
 ./r.stream_power.py dem=elevation time=100 outlets=334736.07613,4747333.21315
