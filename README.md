# MabSim
A monoclonal antibody iIoT simulation using Node-RED and Python

## Install

Install node and node-red.

Clone this repository into ~/.node-red.

## Run

On the command line:
```
cd ~.node-red/mabsim
node-red
```

Automatically will start in http://127.0.0.1:1880

## Security

After installation run: 
```
node-red-admin hash-pw
```
Change authentitation settings at ~/.node-red/settings.js by adding the produced hash. New users can be created with different privileges here aswell.
