#!/bin/bash
set -e

# Namespace
NAMESPACE="omnibioai"

echo "Creating namespace (if not exists)..."
kubectl get namespace $NAMESPACE >/dev/null 2>&1 || kubectl create namespace $NAMESPACE

echo "Applying Kubernetes manifests..."
kubectl apply -f ../kubernetes/ -n $NAMESPACE

echo "Waiting for pods to be ready..."
kubectl wait --for=condition=Ready pod -l app=omnibioai -n $NAMESPACE --timeout=300s
kubectl wait --for=condition=Ready pod -l app=ollama -n $NAMESPACE --timeout=300s
kubectl wait --for=condition=Ready pod -l app=redis -n $NAMESPACE --timeout=300s
kubectl wait --for=condition=Ready pod -l app=postgres -n $NAMESPACE --timeout=300s

echo "Deployment completed!"
echo "OmniBioAI: $(minikube service omnibioai-service -n $NAMESPACE --url)"
echo "Ollama: $(minikube service ollama-service -n $NAMESPACE --url)"
