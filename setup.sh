#!/bin/sh

if [ ${BASH_VERSINFO[0]} -le 2 ]; then
  echo "[ERROR]: You must source this setup script (not run it in a sub-shell). Use like $ source setup.sh"
  exit 1
fi

if ! type add_to_PATH &> /dev/null; then

### Adapted from https://unix.stackexchange.com/questions/4965/keep-duplicates-out-of-path-on-source
function add_to_PATH () {
  for d; do

    d=$(cd -- "$d" && { pwd -P || pwd; }) 2>/dev/null  # canonicalize symbolic links
    if [ -z "$d" ]; then continue; fi  # skip nonexistent directory

    if [ "$d" == "/usr/bin" ] || [ "$d" == "/usr/bin64" ] || [ "$d" == "/usr/local/bin" ] || [ "$d" == "/usr/local/bin64" ]; then
      case ":$PATH:" in
        *":$d:"*) :;;
        *) export PATH=$PATH:$d;;
      esac
    else
      case ":$PATH:" in
        *":$d:"*) :;;
        *) export PATH=$d:$PATH;;
      esac
    fi
  done
}

fi

if ! type add_to_LD_LIBRARY_PATH &> /dev/null; then

function add_to_LD_LIBRARY_PATH () {
  for d; do

    d=$(cd -- "$d" && { pwd -P || pwd; }) 2>/dev/null  # canonicalize symbolic links
    if [ -z "$d" ]; then continue; fi  # skip nonexistent directory

    if [ "$d" == "/usr/lib" ] || [ "$d" == "/usr/lib64" ] || [ "$d" == "/usr/local/lib" ] || [ "$d" == "/usr/local/lib64" ]; then
      case ":$LD_LIBRARY_PATH:" in
        *":$d:"*) :;;
        *) export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$d;;
      esac
    else
      case ":$LD_LIBRARY_PATH:" in
        *":$d:"*) :;;
        *) export LD_LIBRARY_PATH=$d:$LD_LIBRARY_PATH;;
      esac
    fi
  done
}

fi

if ! type is_in_PATH &> /dev/null; then

function is_in_PATH () {
  case ":$PATH:" in
  *:${1}:*) return 0 ;;
  *) return 1 ;;
  esac
}

fi

if [ "1" == "1" ]; then
  if [ -z "${WCSIMDIR}" ]; then
    echo "[ERROR]: WCSIMDIR is not set."
    return 1
  fi
  add_to_LD_LIBRARY_PATH "${WCSIMDIR}/"
fi

if ! which root-config &> /dev/null; then
  echo "[ERROR]: Cannot find root-config in the PATH, is ROOT set up?"
  return 1
fi

#if [ "$(root-config --version)" != "6.24.06" ]; then
#  echo "[WARNING]: opticalAnalysis was built against ROOT version: 6.24.06 living in /user/software/root/root-6.24.06-x86_64-cc7-48, but the currently set up ROOT appears to be version: $(root-config --version)."
#fi


export OPTICALANALYSIS=/user/sjenkins/HyperK/HK-IDLISOAnalysis/build/Linux

add_to_PATH "${OPTICALANALYSIS}/bin"
add_to_LD_LIBRARY_PATH "${OPTICALANALYSIS}/lib"
