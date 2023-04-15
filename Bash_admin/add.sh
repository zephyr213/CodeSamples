#!/bin/sh

for user in `more ulist`
do
echo "$user"
useradd -g csci440 $user
echo "$user@440" | passwd --stdin "$user"
#chage -d 0 $user
mkdir /home/$user/public_html
chmod 711 /home/$user
chown $user:csci440 /home/$user/public_html
chmod 755 /home/$user/public_html
done
