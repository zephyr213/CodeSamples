#if [ "$1" == "-h" ] ; then
#    echo "create user with home at either /data/USER or /data1/USER"
#    echo "sh adduserhpc.sh [0-->/data or 1-->/data1] [USERNAME] [GROUP] 
#    exit 0
#fi


groupname=$3
username=$2
nfsname=nfs"$1"

if [ "$1" == "0" ]
then
	dirname=/data/"$2"
fi


if [ "$1" == "1" ]
then
        dirname=/data1/"$2"
fi


useradd -m -d ${dirname} -g ${groupname} ${username}

userid=$(id -u $username)

echo $userid 
echo $nfsname
echo $groupname

echo "Now creating users on ${nfsname}"

echo $(ssh ${nfsname} "sh /root/adduser.sh ${userid} ${username} ${groupname}")



echo "Now creating password"

echo $username | passwd $username --stdin

echo "Password set as ${username}"

echo "All done!"
