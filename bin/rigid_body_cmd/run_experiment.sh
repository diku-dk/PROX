#!/bin/bash

for i in $(seq 1 10);
do
  cat > cfgs/2016/funnel_dims_offline_$(hostname)_${i}.cfg << EOF
import = cfgs/2016/funnel_dims_offline_$(hostname).cfg
matlab_file = funnel_dims_profile_${i}.m
EOF
./rigid_body_cmd cfgs/2016/funnel_dims_offline_$(hostname)_${i}.cfg
done
