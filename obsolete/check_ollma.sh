Usage Examples
Ollama server
command: ["./wait-for-service.sh", "ollama", "11434", "60", "--", "gunicorn", "omnibioai.wsgi:application", "--bind", "0.0.0.0:8000", "--workers", "3"]

PostgreSQL
./wait-for-service.sh db 5432 30 -- python manage.py migrate

Redis
./wait-for-service.sh redis 6379 -- echo "Redis is ready"