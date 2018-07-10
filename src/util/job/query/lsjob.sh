#!/bin/bash

query_string="real-prop|eval-tsurff|eval-winop"

for p in $(pgrep "$query_string"); do pwdx $p; done

