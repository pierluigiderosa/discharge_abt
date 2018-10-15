#!/usr/bin/env python

############################################################################
#
# MODULE:      corrivazione.py
# AUTHOR(S):   Pierluigi De Rosa
# PURPOSE:     Calcolo della portata di piena con il metodo del ABD tevere
# COPYRIGHT:   (C) 2018 by Pierluigi De Rosa
#              
#
#              This program is free software under the GNU General Public
#              License (>=v3.0) and comes with ABSOLUTELY NO WARRANTY.
#              See the file COPYING that comes with GRASS
#              for details.
#
#############################################################################

#%Module
#%  description: Calcola area contribuente di monte.
#%  keywords: Ortho
#%End

#%option
#% key: dem
#% type: string
#% gisprompt: old,cell,raster
#% description: DTM raster
#% required : yes
#%end

#%option
#% key: time
#% type: integer
#% description: Returning time
#% options: -0-500
#% answer: 0
#%end


#%option
#% key: outlets
#% type: string
#% description: coordinates of outlets
#% required: yes
#%end

import os
import sys,time,math
import numpy as np
import datetime
import grass.script as grass
from grass.pygrass.vector import VectorTopo
from grass.pygrass.vector.geometry import Point
from grass.pygrass.modules.shortcuts import general as g

if "GISBASE" not in os.environ:
	print "You must be in GRASS GIS to run this program."
	sys.exit(1)

def main():
	dem=options['dem']
	TR=options['time'] #TODO Time of concentration
	outlet=options['outlets']
	outlets = outlet.split(',')
	try:
		TR=int(TR)
	except:
		print 'TR is not a number'
		sys.exit()
	grass.use_temp_region()
	#get region in order to estimate the threshold as 1/1000 of total cells
	grass.run_command('g.region',raster=dem)
	
	regione=grass.region()
	thrshold=float(regione['cells'])/300
	#stream and drainage determination
	grass.run_command('r.watershed', elevation=dem, threshold=700, stream='raster_streams', drainage='drainage',overwrite=True,flags='s')
	
	#the radius is little more than the current resolution
	radius=regione['nsres']*1.4
	grass.run_command('r.circle', output='circle', coordinate=outlet, max=radius,overwrite=True) #%(str(outlets[0]),str(outlets[1]))
	#get the distances and take the shortest distance
	distances=grass.read_command('r.distance', map='circle,raster_streams')
	list_dist=distances.split('\n')
	list_dist.remove('')
	list_tuple=[]
	for distance in list_dist:
		dist=distance.split(':')
		my_tupla=dist[0],dist[1],float(dist[2]),dist[3],dist[4],dist[5],dist[6]
		list_tuple.append(my_tupla)
	tuple_orderedByDistance=sorted(list_tuple, key=lambda distanza: distanza[2])
	del(distances,list_tuple,list_dist)
	#calculate the basin and read its statistics
	outlet=tuple_orderedByDistance[0][-2:]
	xoutlet=float(outlet[0])
	youtlet=float(outlet[1])
	grass.run_command('r.water.outlet',input='drainage',output='basin',coordinates=str(xoutlet)+','+str(youtlet) , overwrite=True)
	statistics=grass.read_command('r.univar',map=dem, zones='basin')
	main_stat=statistics.splitlines()[-9:]
	
	
	#order the stream network
	grass.run_command('r.mask',raster='basin')
	grass.run_command('r.stream.order',stream_rast='raster_streams', direction='drainage', elevation=dem,horton='horton',overwrite=True)
	stream_stat=grass.read_command('r.stream.stats', stream_rast='horton', direction='drainage', elevation=dem,flags='o')
	network_statistics=stream_stat.split('\n')
	network_statistics.remove('')
	#get the max order
	network_statistics[-1].split()
	total_length=float(network_statistics[-1].split(',')[2])
	area_basin=float(network_statistics[-1].split(',')[3])
	area_basin_Ha=area_basin*100
	mean_elev=float(main_stat[3].split(':')[-1])
	min_elev=float(main_stat[0].split(':')[-1])
	max_elev=float(main_stat[1].split(':')[-1])
	deltaH=max_elev-min_elev
	average_slope=float(network_statistics[-1].split(',')[4])
	grass.run_command('r.mask',flags='r')
	
	TcGiandotti=(4*np.sqrt(area_basin)+1.5*total_length)/(0.8*np.sqrt(mean_elev-min_elev))
	g.message("TC using Giandotti...")
	grass.info(str(TcGiandotti))
	
	TcKirpich=0.945*(total_length**3./deltaH)**0.385
	g.message("TC using Kirpich...")
	grass.info(str(TcKirpich))
	if area_basin_Ha > 100: #TODO controlla i riferimenti
		corrivazione = TcGiandotti
		grass.info('using giandotti')
	else:
		corrivazione = TcKirpich
		grass.info('using Kirpich')
	if corrivazione < 24:
		aPar='a24@PERMANENT'
		bPar='b24@PERMANENT'
		kPar='k24@PERMANENT'
	else:
		aPar='a15@PERMANENT'
		bPar='b15@PERMANENT'
		kPar='k15@PERMANENT'
	CNmap = 'CN@PERMANENT'
	
	corrivazione=TcGiandotti
	aStat=grass.read_command('r.univar',map=aPar, zones='basin')
	aMain_stat=aStat.splitlines()[12].split(':')[-1]	
	aMain_stat=float(aMain_stat)
	bStat=grass.read_command('r.univar',map=bPar, zones='basin')
	bMain_stat=bStat.splitlines()[12].split(':')[-1]	
	bMain_stat=float(bMain_stat)
	kStat=grass.read_command('r.univar',map=kPar, zones='basin')
	kMain_stat=kStat.splitlines()[12].split(':')[-1]	
	kMain_stat=float(kMain_stat)
	CNstat = grass.read_command('r.univar',map=CNmap, zones='basin')
	CN=CNstat.splitlines()[12].split(':')[-1]
	CN=float(CN)
	
	g.message('area basin in km2: ')
	print area_basin
	print 'mean elev: '
	print mean_elev-min_elev
	print 'delta H:'
	print deltaH
	print 'total reach length: '
	print total_length
	print 'a mean:'
	print aMain_stat
	print '\n b mean: '
	print bMain_stat
	print '\n k mean: '
	print kMain_stat
	print 'CN mean:'
	print CN
	
	f_K_T = 1-kMain_stat*(0.45+0.799*np.log(-np.log(1-1./TR)))
	print 'f(k,T): '
	print f_K_T
	
	h=aMain_stat*corrivazione**bMain_stat*f_K_T
	print '\n h main:'
	print h
	X1 = 100*corrivazione/(0.236+0.062*corrivazione)
	X2 = 0.003*corrivazione+0.0234
	Pa = 100 - area_basin_Ha/(X1+X2*area_basin_Ha)
	Ha = h*Pa/100
	S1 = (1000/CN)-10
	Pn = (Ha-5.08*S1)**2/(Ha+20.32*S1)
	Qc = (1/360.)*Pn*area_basin_Ha/corrivazione
	
	print 'discharge: '
	print Qc
	
	#print table.columns.types()
	#[u'INTEGER', u'TEXT', u'integer', u'double precision']
	
	
	
	'''
	------------------------------
	START CALCULATION OF LOCAL UPSTREAM SLOPE
	------------------------------
	'''
	#offsets for moving windows
	offsets = [d
		   for j in xrange(1,1+1)
		   for i in [j,-j]
		   for d in [(i,0),(0,i),(i,i),(i,-i)]]
	#rename dtm as elevation for future calculation if not exist
	if not VectorTopo('elevation').exist():
		grass.run_command('g.rename',raster="%s,elevation" % dem)
		elev_renamed=True
	
	#define drainage direction
	drainage_incoming = [2,4,3,1,6,8,7,5]
	drainage_outcoming = []
	diag_dist= (regione['nsres']**2+regione['ewres']**2)**0.5
	# [(1, 0), (0, 1), (1, 1), (1, -1), (-1, 0), (0, -1), (-1, -1), (-1, 1), 
	cell_dists = [regione['nsres'], 
					 regione['ewres'],
					 diag_dist,
					 diag_dist,
					 regione['nsres'],
					 regione['ewres'],
					 diag_dist,
					 diag_dist
					]
	# define the calculation term
	terms = ["(drainage[%d,%d] == %d && not(isnull(raster_streams[0,0])) && not(isnull(raster_streams[%d,%d])) )"
				 % ((offsets[j]+tuple([drainage_incoming[j]])+offsets[j]))
				for j in range(len(drainage_incoming))]
	
	   
	   
	 #define the operation expression
	terms_calc = [ "(elevation[%d,%d] - elevation) * %s" 
				% (offsets[j]+(terms[j],) ) for j in range(len(terms))]
	
	terms_calc_slope = [ "( (elevation[%d,%d] - elevation)/%10.4f ) * %s" 
				% (offsets[j]+(cell_dists[j],)+(terms[j],)) for j in range(len(terms))]
	
	expr = "num_cells_drain_into = (%s)" % " + ".join(terms)
	expr1 = "elevation_percentile4 = if(isnull(raster_streams),null(),(%s))" % " + ".join(terms)
	expr2 = "elevdiff_drain_into = %s" % " + ".join(terms_calc)
	expr3 = "slope_drain_into = %s" % " + ".join(terms_calc_slope)       
	
	# do the r.mapcalc calculation with the moving window
	# exclude the num_cell_calculation_into
	#grass.mapcalc( expr )
	#print expr2
	#grass.mapcalc(  expr2 , overwrite=True)
	#print expr3
	#grass.mapcalc(  expr3 , overwrite=True)    
	
	'''
	------------------------------
	START CALCULATION OF 2KM UPSTREAM SLOPE
	------------------------------
	'''
	#create an outlet vector
	new = VectorTopo('outlet')
	COLS = [(u'cat',       'INTEGER PRIMARY KEY')]
	new.open('w', tab_name='outlet', tab_cols=COLS)
	new.write(Point( xoutlet , youtlet ), cat=1, )
	new.table.conn.commit()
	new.table.execute().fetchall()
	new.close()
			

	
	
	new = VectorTopo('output')
	COLS = [(u'cat',       'INTEGER PRIMARY KEY'), (u'discharge',    u'double precision')
		]
		
	if new.exist():
		g.message('The vector exist: it will be named as')
		new_name = 'output' + str(datetime.datetime.now().time())[:8].replace(':','_')
		grass.info( new_name  )
		new = VectorTopo( new_name )
		new.open('w', tab_name=new_name, tab_cols=COLS)
	else:
		new.open('w', tab_name='output', tab_cols=COLS)
		
	
	
	new.write(Point( xoutlet , youtlet ), (float(Qc),) )
	new.table.conn.commit()
	new.table.execute().fetchall()
	new.close()
	
	#cleaning part
	if elev_renamed:
		grass.run_command('g.rename',raster='elevation,%s' % dem)
	grass.del_temp_region()
	grass.run_command('g.remove',flags='f', type='raster', name='main_stream,basin,circle,drainage,horton,raster_streams,slope_drain_into')
	
	grass.run_command('g.remove',flags='f', type='vector', name='main_stream,nodes,outlet')
	grass.run_command('g.remove',type='vector',pattern='main_stream*',flags='f')

	    
	
if __name__ == "__main__":
	options, flags = grass.parser()
	sys.exit(main())
