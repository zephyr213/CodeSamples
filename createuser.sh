# Create new user script
# Getting user basics 

echo Username?
read VarUname
echo Groups?
read VarGname
echo Full Name?
read VarFname

echo $VarUname
echo $VarGname

# generate homepath
homepath=/data1/$VarGname/$VarUname
echo $homepath


# Create user on cluster
useradd -m -d ${homepath} -g ${VarGname} -c '"${VarFname}"' ${VarUname}

# generate default password
echo "$VarUname@hpc" | passwd --stdin "$VarUname"


# getting user id
Vuid=$(id -u ${VarUname})
echo $Vuid

# create user on NFS server with the same user id and group.
sshpass -p "Penguin" ssh nfs1 "useradd -u $Vuid -g $VarGname -s /sbin/nologin $VarUname"
