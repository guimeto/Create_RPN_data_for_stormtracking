r.compile -src libtracks.ftn -includes include -librmn 
ar rv lib/libtracks.a libtracks.o 
ranlib lib/libtracks.a 


#r.build -src make_tracks_copy_uv.ftn -includes include -libpath lib 
#rm -f make_tracks_copy_uv.o  make_tracks_copy_uv.f
#r.compile -libpath ./lib  /HOME/emmanuel/tracking/lib -src make_tracks_copy_uv.ftn -o make_tracks_copy_uv -includes include -libappl "tracks mdp2_cmq  stats 3d div rmn_013"

rm -f make_tracks_narr.o  make_tracks_narr.f
r.build -src make_tracks_narr.ftn -includes include -libpath lib 
r.compile -libpath ./lib  /HOME/emmanuel/tracking/lib -src make_tracks_narr.ftn -o make_tracks_narr -includes include -libappl "tracks mdp2_cmq  stats 3d div rmn_013"
#

