import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

const API = `http://127.0.0.1:${process.env.VENOMGAP_API_PORT ?? 8010}`;

export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    // The API holds the model in memory; proxying keeps the browser on one origin.
    // Port 8000 is a crowded default on a development machine, so the API defaults to 8010.
    // Set VENOMGAP_API_PORT when a supervisor has placed the API somewhere else.
    proxy: { "/api": API, "/health": API },
  },
  build: { outDir: "dist", sourcemap: true },
});
