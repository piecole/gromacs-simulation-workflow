#!/bin/bash

# Create a temporary file for the filter
cat > /tmp/email-filter.sh << 'EOF'
#!/bin/bash
sed '/[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}/d'
EOF

chmod +x /tmp/email-filter.sh

# Run filter-branch
git filter-branch --force --index-filter \
    'git ls-files | while read file; do
        if [[ "$file" =~ \.(py|txt|md|sh)$ ]]; then
            git show ":$file" | /tmp/email-filter.sh > "$file"
            git add "$file"
        fi
    done' \
    --prune-empty --tag-name-filter cat -- --all

# Clean up
rm -rf /tmp/email-filter.sh 