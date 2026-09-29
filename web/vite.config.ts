import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    // The API holds the model in memory; proxying keeps the browser on one origin.
    proxy: { "/api": "http://127.0.0.1:8000", "/health": "http://127.0.0.1:8000" },
  },
  build: { outDir: "dist", sourcemap: true },
});
