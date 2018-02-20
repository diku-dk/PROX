#
# We do a sleep command after each launch to avoid problems with reading/writing TetGen
# files simultaneously. The sleep time should be large enough to make sure that any past
# launches is done setting up their scenes using TetGen.
#
#
# Use the jobs command to see what jobs are still running
#
pwd
echo 
echo 
echo 'Tower is running now'
nohup ../rigid_body_cmdD 2017_SCA_tower.cfg
echo 
echo 
echo 'Wall is running now'
nohup ../rigid_body_cmdD 2017_SCA_wall.cfg
echo 
echo 
echo 'Jacks in Funnel is running now'
nohup ../rigid_body_cmdD 2017_SCA_jacks_in_funnel.cfg
echo 
echo 
echo 'Glasses in glass is running now'
nohup ../rigid_body_cmdD 2017_SCA_glasses_in_glass.cfg
#echo 
#echo 
#echo 'Colosseum is running now'
#nohup ../rigid_body_cmdD 2017_SCA_colosseum.cfg
echo 
echo 
echo 'Jacks in glass is running now'
nohup ../rigid_body_cmdD 2017_SCA_jacks_in_glass.cfg
echo 
echo 
echo 'Spheres in glass is running now'
nohup ../rigid_body_cmdD 2017_SCA_spheres_in_glass.cfg
echo 
echo 
echo 'Shoot is running now'
nohup ../rigid_body_cmdD 2017_SCA_shoot.cfg
echo 
echo 
echo 'Pantheon is running now'
nohup ../rigid_body_cmdD 2017_SCA_pantheon.cfg
#echo 
#echo 
#echo 'Temple is running now'
#nohup ../rigid_body_cmdD 2017_SCA_temple.cfg

