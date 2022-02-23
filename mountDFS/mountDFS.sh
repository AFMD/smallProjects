#!/usr/bin/env bash

# mounts the physics DFS folder, use with sudo, will prompt for

echo Enter your physics username:

read PHYSICS_USER

THIS_USER=`pstree -lu -s $$ | grep --max-count=1 -o '([^)]*)' | head -n 1 | sed 's/[()]//g'`

MOUNT_POINT=/mnt/dfs
mkdir -p "${MOUNT_POINT}"

CIFS_SERVER="dc6.physics.ox.ac.uk"
#CIFS_SERVER="physics.ox.ac.uk"

while : ; do
  umount "${MOUNT_POINT}" > /dev/null 2>&1 || true
  mount -t cifs //${CIFS_SERVER}/dfs /mnt/dfs -o username=${PHYSICS_USER},workgroup=PHYSICS,vers=2.1,uid=`id -u ${THIS_USER}`,iocharset=utf8,gid=`id -g ${THIS_USER}`
  #cd "${MOUNT_POINT}/home/${PHYSICS_USER}" > /dev/null 2>&1
  cd "${MOUNT_POINT}/home/${PHYSICS_USER}"
  if [ $? -eq 0 ]; then
    break
  fi
done
echo "Mounted sucessfully at ${MOUNT_POINT}"
