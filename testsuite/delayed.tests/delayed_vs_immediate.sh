#!/bin/sh
#
# The delayed update must walk the same Markov chain as the immediate one.
#
# Tests 36 to 38 drive delayed_update_mod's own routines. What they cannot
# reach is Upgrade2's delayed branch, which links those routines to the
# Woodbury chain: the reconstructed rows feed v, the reconstructed columns feed
# the rank-d update, and a mistake there produces a plausible chain rather
# than an obviously wrong one.
#
# Here a sampler is run twice on each parameter set, once with the delay off and
# once with it on, from the same seeds. We require the auxiliary field
# configuration left behind to be identical. Up to Metropolis near-ties the
# two schemes accept exactly the same flips, so this is an equality and not a
# tolerance. An HDF5 build writes confout_0.h5 instead of confout_0; two such
# files holding identical data still differ byte-for-byte (HDF5 tracks
# per-object timestamps by default), so that variant is compared
# dataset-by-dataset with h5diff rather than cmp.
#
# Usage: delayed_vs_immediate.sh <ALF.out> <source dir> <work dir>

set -eu

exe=$1
src=$2
work=$3

if [ ! -x "$exe" ]; then
   echo "SKIP: $exe not built"
   exit 77
fi

status=0

for model in hubbard tv; do
   for arm in immediate delayed; do
      dir="$work/$model.$arm"
      rm -rf "$dir"
      mkdir -p "$dir"
      cp "$src/parameters_$model" "$dir/parameters"
      cp "$src/seeds" "$dir/seeds"
      # The depth is the smallest the panels are exercised at. It has to divide
      # the run badly rather than well: with k = 8 and rank-2 vertices the tV
      # arm ends slices holding a partial panel, so the trailing flush in
      # delay_close is on the path too and not just the periodic one.
      if [ "$arm" = delayed ]; then
         ALF_DELAY_K=8 ; export ALF_DELAY_K
      else
         unset ALF_DELAY_K || true
      fi
      ( cd "$dir" && "$exe" > run.log 2>&1 ) || {
         echo "FAIL: $model/$arm did not run; tail of $dir/run.log:"
         tail -20 "$dir/run.log"
         status=1
         continue
      }
   done

   imm="$work/$model.immediate/confout_0"
   del="$work/$model.delayed/confout_0"
   if [ ! -f "$imm" ] && [ -f "$imm.h5" ]; then
      imm="$imm.h5"
      del="$del.h5"
   fi

   if [ ! -f "$imm" ] || [ ! -f "$del" ]; then
      echo "FAIL: $model: no confout_0 to compare"
      status=1
      continue
   fi

   # Guard against a vacuous pass: if the delayed run did not actually take the
   # delayed path, the two runs are the same run and the comparison proves
   # nothing.
   depth=$(awk '/Delay depth  /{print $(NF)}' "$work/$model.delayed/info")
   if [ "${depth:-0}" -ne 8 ]; then
      echo "FAIL: $model: delayed run reports delay depth '${depth:-<none>}', expected 8"
      status=1
      continue
   fi

   diff_log=""
   case "$imm" in
      *.h5)
         diff_log="$work/$model.h5diff.log"
         if h5diff "$imm" "$del" > "$diff_log" 2>&1; then identical=0; else identical=1; fi
         ;;
      *)
         if cmp -s "$imm" "$del"; then identical=0; else identical=1; fi
         ;;
   esac

   if [ "$identical" -eq 0 ]; then
      echo "PASS: $model: confout_0 identical, delay depth $depth"
   else
      echo "FAIL: $model: confout_0 differs between immediate and delayed"
      [ -n "$diff_log" ] && tail -20 "$diff_log"
      status=1
   fi
done

exit $status