#!/bin/bash
   
matlab -batch "neu_sim(prop_formation_rate=${1},separation_distance=${2},sites=${3})"
