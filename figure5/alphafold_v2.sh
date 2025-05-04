###code for alphafold v2

python3 docker/run_docker.py \
    --fasta_paths="$fasta_path" \
    --max_template_date=2020-05-14 \
    --enable_gpu_relax=true \
    --model_preset=multimer \
    --data_dir="$base_output_dir/data" \
    --output_dir="$output_dir" \
    --num_multimer_predictions_per_model=1