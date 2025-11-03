# Deploying to Vercel

This Flask application is configured to run on Vercel as a serverless function.

## Important Limitations

⚠️ **Vercel has several limitations for this application:**

1. **Function Timeout**: 
   - Free tier: 10 seconds
   - Pro tier: 60 seconds (configurable up to 300 seconds in `vercel.json`)
   - Video processing may timeout on longer videos

2. **Request Size Limit**: 
   - Maximum 4.5MB for request body
   - For larger files, consider chunked uploads or external storage

3. **Temporary Storage**:
   - Files are stored in `/tmp` which is ephemeral
   - Files are deleted after function execution completes
   - Maximum 512MB in `/tmp`

4. **Model Weights**:
   - Large ML model weights must be included in deployment
   - This may increase deployment size significantly

## Deployment Steps

1. **Install Vercel CLI** (if not already installed):
   ```bash
   npm i -g vercel
   ```

2. **Login to Vercel**:
   ```bash
   vercel login
   ```

3. **Deploy**:
   ```bash
   vercel
   ```

4. **For Production**:
   ```bash
   vercel --prod
   ```

## Configuration Files

- `vercel.json`: Vercel configuration with Python runtime and extended timeout
- `.vercelignore`: Files to exclude from deployment
- `pyproject.toml`: Fixed to resolve setuptools errors

## Alternative Deployment Options

If you encounter timeout or size issues, consider:

1. **Railway**: Better for long-running processes
2. **Render**: Supports Flask with persistent storage
3. **Fly.io**: Good for containerized apps with large dependencies
4. **AWS Lambda + API Gateway**: Serverless with longer timeouts (up to 15 minutes)

## Troubleshooting

### Build Errors

If you see setuptools errors, ensure `pyproject.toml` is properly configured (already fixed).

### Timeout Errors

- Reduce video duration for testing
- Process videos in chunks
- Use Vercel Pro plan for longer timeouts

### Size Errors

- Model weights may be too large
- Consider using Vercel's large file handling or external storage

