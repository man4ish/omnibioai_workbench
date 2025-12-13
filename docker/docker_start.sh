RUN service mysql start && \
    mysql -u root -proot -e "USE omnibioai; \
    CREATE TABLE IF NOT EXISTS workflow_runs ( \
        id INT AUTO_INCREMENT PRIMARY KEY, \
        workflow_name VARCHAR(255) NOT NULL, \
        engine VARCHAR(50) NOT NULL, \
        entrypoint TEXT NOT NULL, \
        parameters JSON, \
        start_time DATETIME NOT NULL, \
        end_time DATETIME DEFAULT NULL, \
        state VARCHAR(20) NOT NULL, \
        outputs JSON, \
        logs LONGTEXT \
    );"